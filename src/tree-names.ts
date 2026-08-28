/**
 * Produce Newick-safe, stable and unique tip identifiers while preserving row
 * order. The numeric FastTree input is restored through this exact mapping in
 * both browser and CLI builds.
 */
export function treeTipNames(names: readonly string[]): string[] {
  const used = new Set<string>();
  return names.map((name, index) => {
    const base = name.replace(/[^A-Za-z0-9_.|*+\-]/g, "_") || `tip_${index + 1}`;
    let candidate = base, suffix = 1;
    while (used.has(candidate)) candidate = `${base}__${index + 1}_${suffix++}`;
    used.add(candidate); return candidate;
  });
}
