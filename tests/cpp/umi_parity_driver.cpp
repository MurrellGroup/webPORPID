#include "webporpid/binary.hpp"
#include "webporpid/core.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

int main() {
  const std::vector<std::pair<std::string, std::uint32_t>> counts{
    {"AACCGGT", 3}, {"AACCGGTA", 4}, {"AACCGGTT", 100}, {"AACCGTT", 2},
    {"AAACCGGTT", 2}, {"TTGGCCAA", 20}, {"TTGGCCAT", 1},
  };
  std::vector<std::uint8_t> encoded; webporpid::binary::magic(encoded, "WPN1");
  webporpid::binary::number(encoded, static_cast<std::uint32_t>(counts.size()));
  for (const auto& [tag, count] : counts) {
    webporpid::binary::number(encoded, std::uint16_t{0}); webporpid::binary::string(encoded, tag);
    webporpid::binary::number(encoded, count);
  }
  webporpid::Config config; config.samples.push_back(webporpid::Sample{});
  config.parameters.lda_threshold = 0.995; config.parameters.family_size_threshold = 1;
  std::string error; const auto model = webporpid::build_family_model(encoded, config, error);
  if (!error.empty()) { std::cerr << error << '\n'; return 1; }
  std::cout << std::setprecision(17);
  for (const auto& row : model) std::cout << row.umi << '\t' << row.parent << '\t' << row.family_size << '\t'
    << row.probability << '\t' << static_cast<int>(row.disposition) << '\n';
}
