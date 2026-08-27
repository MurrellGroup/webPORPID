#include "webporpid/core.hpp"
#include "webporpid/binary.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>

namespace webporpid {
namespace {

using Counts = std::map<std::uint16_t, std::map<std::string, std::uint32_t>>;
using SparseRow = std::vector<std::pair<std::size_t, double>>;

bool decode_cutoffs(std::span<const std::uint8_t> bytes, std::vector<std::uint64_t>& output, std::string& error) {
  binary::Reader reader(bytes); std::uint32_t count = 0;
  if (!reader.magic("WPT1") || !reader.number(count) || count > 65535) {
    error = "The downsampling threshold table is invalid."; return false;
  }
  output.resize(count);
  for (auto& threshold : output) if (!reader.number(threshold)) { error = "The threshold table is truncated."; return false; }
  if (!reader.done()) { error = "The threshold table has trailing bytes."; return false; }
  return true;
}

bool selected(const SpoolRecord& record, const std::vector<std::uint64_t>& cutoffs) {
  return record.sample >= cutoffs.size() || record.sampling_hash <= cutoffs[record.sample];
}

bool decode_counts(std::span<const std::uint8_t> bytes, Counts& output, std::string& error) {
  binary::Reader reader(bytes); std::uint32_t count = 0;
  if (!reader.magic("WPN1") || !reader.number(count) || count > 100000000) {
    error = "The merged UMI count table is invalid."; return false;
  }
  for (std::uint32_t index = 0; index < count; ++index) {
    std::uint16_t sample = 0; std::string umi; std::uint32_t family_size = 0;
    if (!reader.number(sample) || !reader.string(umi) || !reader.number(family_size)) {
      error = "The merged UMI count table is truncated."; return false;
    }
    output[sample][umi] += family_size;
  }
  if (!reader.done()) { error = "The merged UMI count table has trailing bytes."; return false; }
  return true;
}

std::vector<std::pair<std::string, double>> one_edit_words(std::string_view tag, bool include_identity_mutations = false) {
  constexpr double error = 0.005;
  constexpr double insertion = error * 0.4 * 0.25;
  constexpr double deletion = error * 0.4;
  constexpr double mutation = error * 0.2 / 3.0;
  std::vector<std::pair<std::string, double>> output;
  output.reserve(4 * (tag.size() + 1) + tag.size() + 3 * tag.size());
  for (char base : std::string_view("ACTG")) for (std::size_t position = 0; position <= tag.size(); ++position) {
    std::string word(tag); word.insert(word.begin() + static_cast<std::ptrdiff_t>(position), base);
    output.emplace_back(std::move(word), insertion);
  }
  for (std::size_t position = 0; position < tag.size(); ++position) {
    std::string word(tag); word.erase(position, 1); output.emplace_back(std::move(word), deletion);
  }
  for (char base : std::string_view("ACTG")) for (std::size_t position = 0; position < tag.size(); ++position) {
    if (!include_identity_mutations && tag[position] == base) continue;
    std::string word(tag); word[position] = base; output.emplace_back(std::move(word), mutation);
  }
  return output;
}

std::map<std::size_t, double> one_step(std::string_view observed,
    const std::unordered_map<std::string, std::size_t>& candidates) {
  std::map<std::size_t, double> likelihood;
  for (const auto& [word, weight] : one_edit_words(observed)) {
    if (const auto found = candidates.find(word); found != candidates.end()) likelihood[found->second] += weight;
  }
  if (const auto found = candidates.find(std::string(observed)); found != candidates.end())
    likelihood[found->second] = std::pow(1.0 - 0.005, static_cast<double>(observed.size()));
  return likelihood;
}

SparseRow two_edit_likelihood(std::string_view observed,
    const std::unordered_map<std::string, std::size_t>& candidates) {
  std::map<std::size_t, double> likelihood;
  std::unordered_map<std::string, std::map<std::size_t, double>> memo;
  for (const auto& [intermediate, first_weight] : one_edit_words(observed, true)) {
    auto [iterator, inserted] = memo.try_emplace(intermediate);
    if (inserted) iterator->second = one_step(intermediate, candidates);
    for (const auto& [index, probability] : iterator->second) likelihood[index] += first_weight * probability;
  }
  if (const auto found = candidates.find(std::string(observed)); found != candidates.end())
    likelihood[found->second] = std::pow(1.0 - 0.005, static_cast<double>(observed.size()));
  return {likelihood.begin(), likelihood.end()};
}

std::vector<std::pair<std::size_t, double>> lda(const std::vector<SparseRow>& likelihoods,
                                                const std::vector<std::uint32_t>& counts) {
  constexpr double concentration = 0.5, epsilon = 1e-17;
  const std::size_t count = likelihoods.size();
  std::vector<double> prior(count, 1.0 / static_cast<double>(count)), posterior(count, 0.0);
  for (std::size_t iteration = 0; iteration < 1000; ++iteration) {
    for (std::size_t observed = 0; observed < count; ++observed) {
      double total = 0.0;
      for (const auto& [real, probability] : likelihoods[observed]) total += prior[real] * probability;
      if (total == 0.0) continue;
      const double normalized = static_cast<double>(counts[observed]) / total;
      for (const auto& [real, probability] : likelihoods[observed]) posterior[real] += prior[real] * probability * normalized;
    }
    double total = 0.0; for (double value : posterior) total += value + concentration;
    bool converged = true;
    for (std::size_t real = 0; real < count; ++real) {
      const double next = (posterior[real] + concentration) / total;
      converged = converged && std::abs(next - prior[real]) < epsilon / static_cast<double>(count);
      prior[real] = next; posterior[real] = 0.0;
    }
    if (converged) break;
  }
  std::vector<std::pair<std::size_t, double>> output; output.reserve(count);
  for (std::size_t observed = 0; observed < count; ++observed) {
    double best = 0.0, total = 0.0; std::size_t best_index = observed;
    for (const auto& [real, probability] : likelihoods[observed]) {
      const double value = prior[real] * probability; total += value;
      if (value > best) { best = value; best_index = real; }
    }
    output.emplace_back(best_index, total > 0.0 ? best / total : 0.0);
  }
  return output;
}

} // namespace

std::vector<std::uint8_t> partition_sample_hashes(std::span<const std::uint8_t> bytes,
                                                  const Config& config, std::string& error) {
  const auto records = decode_spool(bytes, error); if (!error.empty()) return {};
  std::vector<std::uint64_t> counts(config.samples.size(), 0);
  for (const auto& record : records) if (record.sample < counts.size()) ++counts[record.sample];
  std::vector<std::uint8_t> output; binary::magic(output, "WPS1");
  binary::number(output, static_cast<std::uint32_t>(counts.size()));
  for (auto count : counts) binary::number(output, count);
  return output;
}

std::vector<std::uint8_t> count_families(std::span<const std::uint8_t> bytes,
                                         std::span<const std::uint8_t> cutoff_bytes,
                                         const Config&, std::string& error) {
  std::vector<std::uint64_t> cutoffs;
  if (!decode_cutoffs(cutoff_bytes, cutoffs, error)) return {};
  const auto records = decode_spool(bytes, error); if (!error.empty()) return {};
  Counts counts;
  for (const auto& record : records) if (selected(record, cutoffs)) ++counts[record.sample][record.umi];
  std::size_t entries = 0; for (const auto& [_, families] : counts) entries += families.size();
  std::vector<std::uint8_t> output; binary::magic(output, "WPN1");
  binary::number(output, static_cast<std::uint32_t>(entries));
  for (const auto& [sample, families] : counts) for (const auto& [umi, count] : families) {
    binary::number(output, sample); binary::string(output, umi); binary::number(output, count);
  }
  return output;
}

std::vector<FamilyDecision> build_family_model(std::span<const std::uint8_t> merged_counts,
                                               const Config& config, std::string& error) {
  Counts by_sample; if (!decode_counts(merged_counts, by_sample, error)) return {};
  std::vector<FamilyDecision> output;
  for (const auto& [sample_index, sample_counts] : by_sample) {
    if (sample_index >= config.samples.size() || sample_counts.empty()) continue;
    std::vector<std::string> tags; std::vector<std::uint32_t> counts;
    tags.reserve(sample_counts.size()); counts.reserve(sample_counts.size());
    for (const auto& [tag, count] : sample_counts) {
      if (tag == "REJECTS") {
        output.push_back(FamilyDecision{sample_index, tag, "REJECTED", count,
          1.0 - std::exp(-5.123456789), FamilyDisposition::bpb_reject});
      } else { tags.push_back(tag); counts.push_back(count); }
    }
    if (tags.empty()) continue;
    std::unordered_map<std::string, std::size_t> indices;
    for (std::size_t index = 0; index < tags.size(); ++index) indices.emplace(tags[index], index);
    std::vector<SparseRow> likelihoods; likelihoods.reserve(tags.size());
    for (const auto& tag : tags) likelihoods.push_back(two_edit_likelihood(tag, indices));
    const auto assignments = lda(likelihoods, counts);
    const auto override = config.samples[sample_index].family_size_override;
    const auto family_threshold = override >= 0 ? static_cast<std::uint32_t>(override) : config.parameters.family_size_threshold;
    for (std::size_t observed = 0; observed < tags.size(); ++observed) {
      FamilyDisposition disposition = FamilyDisposition::likely_real;
      if (assignments[observed].second < config.parameters.lda_threshold) disposition = FamilyDisposition::lda_reject;
      else if (tags[observed].size() != 8) disposition = FamilyDisposition::umi_length_reject;
      else if (counts[observed] < family_threshold) disposition = FamilyDisposition::family_size_reject;
      output.push_back(FamilyDecision{sample_index, tags[observed], tags[assignments[observed].first], counts[observed],
                                      assignments[observed].second, disposition});
    }
  }
  return output;
}

std::vector<std::uint8_t> encode_family_model(const std::vector<FamilyDecision>& model) {
  std::vector<std::uint8_t> output; binary::magic(output, "WPM1");
  binary::number(output, static_cast<std::uint32_t>(model.size()));
  for (const auto& decision : model) {
    binary::number(output, decision.sample); binary::string(output, decision.umi); binary::string(output, decision.parent);
    binary::number(output, decision.family_size); binary::number(output, decision.probability);
    binary::number(output, static_cast<std::uint8_t>(decision.disposition));
  }
  return output;
}

std::vector<FamilyDecision> decode_family_model(std::span<const std::uint8_t> bytes, std::string& error) {
  binary::Reader reader(bytes); std::uint32_t count = 0;
  if (!reader.magic("WPM1") || !reader.number(count) || count > 100000000) {
    error = "The family model is invalid."; return {};
  }
  std::vector<FamilyDecision> output; output.reserve(count);
  for (std::uint32_t index = 0; index < count; ++index) {
    FamilyDecision value; std::uint8_t disposition = 0;
    if (!reader.number(value.sample) || !reader.string(value.umi) || !reader.string(value.parent)
        || !reader.number(value.family_size) || !reader.number(value.probability) || !reader.number(disposition)
        || disposition > static_cast<std::uint8_t>(FamilyDisposition::family_size_reject)) {
      error = "The family model is truncated."; return {};
    }
    value.disposition = static_cast<FamilyDisposition>(disposition); output.push_back(std::move(value));
  }
  if (!reader.done()) { error = "The family model has trailing bytes."; return {}; }
  return output;
}

} // namespace webporpid
