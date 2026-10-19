#include "kg/Expression.h"

#include <cctype>
#include <cmath>
#include <map>
#include <stdexcept>

namespace kg {

namespace {

struct Lex {
  enum Kind { Num, Ident, Op, LParen, RParen, Comma, End } kind;
  double num = 0;
  std::string str;
};

std::vector<Lex> lex(const std::string &s) {
  std::vector<Lex> out;
  size_t i = 0;
  while (i < s.size()) {
    const char c = s[i];
    if (std::isspace(static_cast<unsigned char>(c))) {
      ++i;
    } else if (std::isdigit(static_cast<unsigned char>(c)) || (c == '.' && i + 1 < s.size() && std::isdigit(static_cast<unsigned char>(s[i + 1])))) {
      size_t len = 0;
      const double v = std::stod(s.substr(i), &len);
      out.push_back({Lex::Num, v, ""});
      i += len;
    } else if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
      size_t j = i;
      while (j < s.size() && (std::isalnum(static_cast<unsigned char>(s[j])) || s[j] == '_' || s[j] == '.')) ++j;
      out.push_back({Lex::Ident, 0, s.substr(i, j - i)});
      i = j;
    } else if (c == '(') {
      out.push_back({Lex::LParen, 0, "("});
      ++i;
    } else if (c == ')') {
      out.push_back({Lex::RParen, 0, ")"});
      ++i;
    } else if (c == ',') {
      out.push_back({Lex::Comma, 0, ","});
      ++i;
    } else if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^') {
      out.push_back({Lex::Op, 0, std::string(1, c)});
      ++i;
    } else {
      throw std::invalid_argument("expression '" + s + "': unexpected character '" + std::string(1, c) + "'");
    }
  }
  out.push_back({Lex::End, 0, ""});
  return out;
}

} // namespace

// Recursive-descent parser producing RPN.
Expression::Expression(const std::string &text) : text_(text) {
  const std::vector<Lex> toks = lex(text);
  size_t pos = 0;
  static const std::map<std::string, std::pair<Op, int>> funcs = {
      {"min", {Op::Min, 2}}, {"max", {Op::Max, 2}}, {"sqrt", {Op::Sqrt, 1}}, {"abs", {Op::Abs, 1}},
      {"log", {Op::Log, 1}}, {"exp", {Op::Exp, 1}}, {"pow", {Op::PowF, 2}}};
  auto fail = [&](const std::string &msg) { throw std::invalid_argument("expression '" + text + "': " + msg); };

  std::function<void()> parseExpr, parseTerm, parseFactor, parseUnary, parsePrimary;
  parsePrimary = [&]() {
    const Lex &t = toks[pos];
    if (t.kind == Lex::Num) {
      rpn_.push_back({Op::Num, t.num, "", -1});
      ++pos;
    } else if (t.kind == Lex::Ident) {
      ++pos;
      if (toks[pos].kind == Lex::LParen) {
        auto f = funcs.find(t.str);
        if (f == funcs.end()) fail("unknown function '" + t.str + "'");
        ++pos;
        int nargs = 0;
        if (toks[pos].kind != Lex::RParen) {
          while (true) {
            parseExpr();
            ++nargs;
            if (toks[pos].kind == Lex::Comma) {
              ++pos;
              continue;
            }
            break;
          }
        }
        if (toks[pos].kind != Lex::RParen) fail("expected ')'");
        ++pos;
        if (nargs != f->second.second) fail("function '" + t.str + "' takes " + std::to_string(f->second.second) + " argument(s)");
        rpn_.push_back({f->second.first, 0, "", -1});
      } else {
        rpn_.push_back({Op::Var, 0, t.str, -1});
      }
    } else if (t.kind == Lex::LParen) {
      ++pos;
      parseExpr();
      if (toks[pos].kind != Lex::RParen) fail("expected ')'");
      ++pos;
    } else {
      fail("unexpected token '" + t.str + "'");
    }
  };
  parseUnary = [&]() {
    if (toks[pos].kind == Lex::Op && (toks[pos].str == "-" || toks[pos].str == "+")) {
      const bool neg = toks[pos].str == "-";
      ++pos;
      parseUnary();
      if (neg) rpn_.push_back({Op::Neg, 0, "", -1});
    } else {
      parsePrimary();
      if (toks[pos].kind == Lex::Op && toks[pos].str == "^") { // right associative
        ++pos;
        parseUnary();
        rpn_.push_back({Op::Pow, 0, "", -1});
      }
    }
  };
  parseFactor = parseUnary;
  parseTerm = [&]() {
    parseFactor();
    while (toks[pos].kind == Lex::Op && (toks[pos].str == "*" || toks[pos].str == "/")) {
      const Op op = toks[pos].str == "*" ? Op::Mul : Op::Div;
      ++pos;
      parseFactor();
      rpn_.push_back({op, 0, "", -1});
    }
  };
  parseExpr = [&]() {
    parseTerm();
    while (toks[pos].kind == Lex::Op && (toks[pos].str == "+" || toks[pos].str == "-")) {
      const Op op = toks[pos].str == "+" ? Op::Add : Op::Sub;
      ++pos;
      parseTerm();
      rpn_.push_back({op, 0, "", -1});
    }
  };
  parseExpr();
  if (toks[pos].kind != Lex::End) fail("trailing input");
}

std::vector<std::string> Expression::variables() const {
  std::vector<std::string> out;
  for (const auto &t : rpn_)
    if (t.op == Op::Var) {
      bool seen = false;
      for (const auto &o : out) seen |= (o == t.name);
      if (!seen) out.push_back(t.name);
    }
  return out;
}

void Expression::bind(const std::function<int(const std::string &)> &resolver) {
  for (auto &t : rpn_)
    if (t.op == Op::Var) {
      t.slot = resolver(t.name);
      if (t.slot < 0) throw std::invalid_argument("expression '" + text_ + "': unknown parameter '" + t.name + "'");
    }
}

double Expression::eval(const std::vector<double> &v) const {
  double stack[64];
  int sp = 0;
  for (const auto &t : rpn_) {
    switch (t.op) {
      case Op::Num: stack[sp++] = t.num; break;
      case Op::Var: stack[sp++] = v.at(t.slot); break;
      case Op::Neg: stack[sp - 1] = -stack[sp - 1]; break;
      case Op::Sqrt: stack[sp - 1] = std::sqrt(stack[sp - 1]); break;
      case Op::Abs: stack[sp - 1] = std::fabs(stack[sp - 1]); break;
      case Op::Log: stack[sp - 1] = std::log(stack[sp - 1]); break;
      case Op::Exp: stack[sp - 1] = std::exp(stack[sp - 1]); break;
      default: {
        const double b = stack[--sp];
        double &a = stack[sp - 1];
        switch (t.op) {
          case Op::Add: a += b; break;
          case Op::Sub: a -= b; break;
          case Op::Mul: a *= b; break;
          case Op::Div: a /= b; break;
          case Op::Pow:
          case Op::PowF: a = std::pow(a, b); break;
          case Op::Min: a = std::min(a, b); break;
          case Op::Max: a = std::max(a, b); break;
          default: break;
        }
      }
    }
    if (sp >= 63) throw std::runtime_error("expression '" + text_ + "' too deeply nested");
  }
  return stack[0];
}

} // namespace kg
