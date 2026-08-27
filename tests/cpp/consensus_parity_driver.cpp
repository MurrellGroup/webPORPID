#include "webporpid/core.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Generator {
  std::uint32_t state = 0x6d2b79f5U;
  std::uint32_t next() { state = state * 1664525U + 1013904223U; return state; }
  std::size_t below(std::size_t limit) { return (next() >> 8) % limit; }
  char base() { return "ACGT"[below(4)]; }
};

std::string template_sequence(Generator& generator, std::size_t length) {
  std::string output; output.reserve(length);
  for (std::size_t index = 0; index < length; ++index) output.push_back(generator.base());
  return output;
}

std::string noisy(Generator& generator, const std::string& source,
                  std::uint32_t deletion, std::uint32_t insertion, std::uint32_t substitution) {
  std::string output; output.reserve(source.size() + source.size() / 20);
  for (char value : source) {
    const auto event = generator.below(10000);
    if (event < deletion) continue;
    if (event < deletion + insertion) output.push_back(generator.base());
    if (event < deletion + insertion + substitution) {
      char replacement = generator.base(); while (replacement == value) replacement = generator.base(); output.push_back(replacement);
    } else output.push_back(value);
  }
  return output;
}

void emit(std::string_view name, const std::vector<std::string>& reads) {
  double minimum = 0.0; std::vector<webporpid::LowAgreementSite> sites;
  const auto consensus = webporpid::family_consensus(reads, minimum, sites);
  std::cout << name << '\t' << minimum << '\t' << consensus << '\n';
}

} // namespace

int main() {
  Generator generator;
  const auto source = template_sequence(generator, 720);
  emit("identical", std::vector<std::string>(5, source));

  std::vector<std::string> substitutions;
  for (int index = 0; index < 9; ++index) substitutions.push_back(noisy(generator, source, 0, 0, 220));
  emit("substitution", substitutions);

  std::vector<std::string> mixed;
  for (int index = 0; index < 11; ++index) mixed.push_back(noisy(generator, source, 120, 120, 280));
  emit("mixed_indel", mixed);

  std::vector<std::string> high;
  for (int index = 0; index < 13; ++index) high.push_back(noisy(generator, source, 350, 350, 400));
  emit("high_indel", high);

  std::vector<std::string> terminal;
  for (int index = 0; index < 9; ++index) {
    auto read = noisy(generator, source, 100, 100, 180);
    if (index < 5) read = std::string("ACG") + read;
    if (index % 3 == 0 && read.size() > 8) read.erase(0, 5);
    terminal.push_back(std::move(read));
  }
  emit("terminal_overhang", terminal);
}
