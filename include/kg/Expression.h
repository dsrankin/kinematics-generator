// Tiny arithmetic expression language used for dependent parameters,
// e.g. "0.5 * X.mass" or "max(mY + 10, 2 * mZ)".
//
// Grammar: numbers, identifiers ([A-Za-z_][A-Za-z0-9_.]*), + - * / ^,
// unary minus, parentheses, and the functions
// min(a,b) max(a,b) sqrt abs log exp pow(a,b).
#ifndef KG_EXPRESSION_H
#define KG_EXPRESSION_H

#include <functional>
#include <string>
#include <vector>

namespace kg {

class Expression {
 public:
  Expression() = default;
  explicit Expression(const std::string &text);

  const std::string &text() const { return text_; }
  bool empty() const { return rpn_.empty(); }

  // Identifiers referenced by the expression (unique, in order of appearance).
  std::vector<std::string> variables() const;

  // Bind identifiers to slots in a value array. Throws if the resolver returns < 0.
  void bind(const std::function<int(const std::string &)> &resolver);

  // Evaluate with bound variable slots.
  double eval(const std::vector<double> &values) const;

 private:
  enum class Op { Num, Var, Add, Sub, Mul, Div, Pow, Neg, Min, Max, Sqrt, Abs, Log, Exp, PowF };
  struct Tok {
    Op op;
    double num = 0;
    std::string name;
    int slot = -1;
  };
  std::string text_;
  std::vector<Tok> rpn_;
};

} // namespace kg

#endif
