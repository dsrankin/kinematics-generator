// ROOT output (optional library kinegen_root).
//
// Trees written:
//   events        one entry per observed event: ids, label, variable-length
//                 object arrays (obj_*), flattened features (obj_feat, n_feat
//                 per object), MET and, optionally, the full truth record (truth_*).
//   tasks         one entry per task: ids, seeds, topology, feature schema,
//                 object-type vocabulary, parameter names/values and the full
//                 task YAML (enough to regenerate every event).
//   environments  one entry per detector environment: id, seed, YAML.
// With legacy_format the original test.cc trees ("test", "info") are written instead.
#ifndef KG_ROOTWRITER_H
#define KG_ROOTWRITER_H

#include <memory>
#include <string>

#include "kg/Event.h"
#include "kg/Task.h"

namespace kg {

struct RootWriterOptions {
  bool write_truth = true;
  bool legacy_format = false;
};

class RootWriter {
 public:
  explicit RootWriter(const std::string &path, RootWriterOptions opt = {});
  ~RootWriter();
  RootWriter(const RootWriter &) = delete;
  RootWriter &operator=(const RootWriter &) = delete;

  void writeTask(const TaskConfig &task);
  void writeEnvironment(const Environment &env);
  void fill(const TaskConfig &task, const TruthEvent &truth, const ObservedEvent &reco);
  void close();

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace kg

#endif
