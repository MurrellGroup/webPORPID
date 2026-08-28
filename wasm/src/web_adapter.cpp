#include "webporpid/core.hpp"

#include <cstdlib>
#include <sstream>

namespace {
webporpid::Config config;
webporpid::Stats stats;
std::vector<webporpid::FamilyDecision> model;
std::vector<std::uint8_t> result;
std::string error;

std::span<const std::uint8_t> input(const std::uint8_t* pointer, std::uint32_t length) {
  return pointer || length == 0 ? std::span(pointer, length) : std::span<const std::uint8_t>{};
}

std::int32_t fail(std::string message) { error = std::move(message); result.clear(); return -1; }

template <class Callable> std::int32_t invoke(Callable callable) {
  error.clear(); result = callable();
  return error.empty() ? static_cast<std::int32_t>(result.size()) : -1;
}
} // namespace

#if defined(__wasm__)
#define WPP_EXPORT __attribute__((visibility("default"))) __attribute__((used))
#else
#define WPP_EXPORT
#endif

extern "C" {

WPP_EXPORT std::uint8_t* wpp_alloc(std::uint32_t size) {
  return static_cast<std::uint8_t*>(std::malloc(size));
}

WPP_EXPORT void wpp_free(void* pointer) { std::free(pointer); }

WPP_EXPORT const char* wpp_version() { return "0.2.0"; }

WPP_EXPORT std::int32_t wpp_init_config(const std::uint8_t* pointer, std::uint32_t length) {
  error.clear(); result.clear(); model.clear(); stats = {};
  if (!pointer && length) return fail("The configuration pointer is null.");
  if (!webporpid::decode_config(input(pointer, length), config, error)) return -1;
  stats.per_sample.assign(config.samples.size(), 0); return static_cast<std::int32_t>(config.samples.size());
}

WPP_EXPORT std::int32_t wpp_preprocess(const std::uint8_t* pointer, std::uint32_t length, std::uint32_t first_ordinal) {
  if (!pointer && length) return fail("The FASTQ batch pointer is null.");
  return invoke([&] { return webporpid::preprocess_batch(config, stats,
    std::string_view(reinterpret_cast<const char*>(pointer), length), first_ordinal); });
}

WPP_EXPORT std::int32_t wpp_partition_counts(const std::uint8_t* pointer, std::uint32_t length) {
  if (!pointer && length) return fail("The spool partition pointer is null.");
  return invoke([&] { return webporpid::partition_sample_hashes(input(pointer, length), config, error); });
}

WPP_EXPORT std::int32_t wpp_count_families(const std::uint8_t* pointer, std::uint32_t length,
                                           const std::uint8_t* thresholds, std::uint32_t thresholds_length) {
  if ((!pointer && length) || (!thresholds && thresholds_length)) return fail("A family-count pointer is null.");
  return invoke([&] { return webporpid::count_families(input(pointer, length), input(thresholds, thresholds_length), config, error); });
}

WPP_EXPORT std::int32_t wpp_build_family_model(const std::uint8_t* pointer, std::uint32_t length) {
  if (!pointer && length) return fail("The merged count pointer is null.");
  error.clear(); model = webporpid::build_family_model(input(pointer, length), config, error);
  if (!error.empty()) { result.clear(); return -1; }
  result = webporpid::encode_family_model(model); return static_cast<std::int32_t>(result.size());
}

WPP_EXPORT std::int32_t wpp_init_family_model(const std::uint8_t* pointer, std::uint32_t length) {
  if (!pointer && length) return fail("The family-model pointer is null.");
  error.clear(); model = webporpid::decode_family_model(input(pointer, length), error);
  return error.empty() ? static_cast<std::int32_t>(model.size()) : -1;
}

WPP_EXPORT std::int32_t wpp_consensus_partition(const std::uint8_t* pointer, std::uint32_t length,
                                                const std::uint8_t* thresholds, std::uint32_t thresholds_length) {
  if ((!pointer && length) || (!thresholds && thresholds_length)) return fail("A consensus pointer is null.");
  return invoke([&] { return webporpid::process_consensus_partition(input(pointer, length),
      input(thresholds, thresholds_length), model, config, error); });
}

WPP_EXPORT std::int32_t wpp_stats() {
  std::ostringstream stream;
  stream << "{\"totalReads\":" << stats.total_reads << ",\"qualityReads\":" << stats.quality_reads
    << ",\"badReads\":" << stats.bad_reads << ",\"shortReads\":" << stats.short_reads
    << ",\"longReads\":" << stats.long_reads << ",\"primerRejects\":" << stats.primer_rejects
    << ",\"idRejects\":" << stats.id_rejects << ",\"demultiplexedReads\":" << stats.demultiplexed_reads
    << ",\"bpbRejects\":" << stats.bpb_rejects << ",\"malformedRecords\":" << stats.malformed_records
    << ",\"perSample\":[";
  for (std::size_t index = 0; index < stats.per_sample.size(); ++index) {
    if (index) stream << ',';
    stream << stats.per_sample[index];
  }
  stream << "]}";
  const auto text = stream.str(); result.assign(text.begin(), text.end()); return static_cast<std::int32_t>(result.size());
}

WPP_EXPORT const std::uint8_t* wpp_result_ptr() { return result.data(); }
WPP_EXPORT std::uint32_t wpp_result_len() { return static_cast<std::uint32_t>(result.size()); }
WPP_EXPORT const char* wpp_error_ptr() { return error.data(); }
WPP_EXPORT std::uint32_t wpp_error_len() { return static_cast<std::uint32_t>(error.size()); }

} // extern "C"
