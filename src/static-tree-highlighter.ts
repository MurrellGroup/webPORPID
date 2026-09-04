import { effectiveAlignment, inspectAlignment, type AlignmentVariant } from "./alignment-utils.ts";
import { functionalSequenceName, uncollapsedSequenceName } from "./sequence-names.ts";
import { ladderizeTree, layoutTree, parseNewick, rootOnOutgroup, type TreeNode } from "./tree.ts";
import { treeTipNames } from "./tree-names.ts";
import type { ResultBundle } from "./types.ts";

const WIDTH = 1600;
const TREE_PANE_WIDTH = 780;
const TREE_BRANCH_WIDTH = 610;
const MATRIX_LEFT = 820;
const MATRIX_WIDTH = 750;
const TOP = 68;
const ROW_HEIGHT = 18;
const BUBBLE_AREA_PER_FAMILY = 25;
const COLORS: Record<string, string> = {
  A: "#f8766d",
  G: "#ffcc00",
  T: "#80e586",
  C: "#a6d5e3",
  U: "#a6d5e3",
  "-": "#9e9e9e",
  N: "#737373",
};

const xml = (value: unknown) => String(value).replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll(">", "&gt;").replaceAll('"', "&quot;");

function leaves(root: TreeNode): string[] {
  const output: string[] = [];
  const visit = (node: TreeNode) => node.children.length ? node.children.forEach(visit) : output.push(node.name);
  visit(root);
  return output;
}

function niceStep(maximum: number): number {
  if (!(maximum > 0)) return 1;
  const target = maximum / 4, power = 10 ** Math.floor(Math.log10(target)), scaled = target / power;
  return (scaled >= 5 ? 5 : scaled >= 2 ? 2 : 1) * power;
}

function familyCounts(bundle: ResultBundle, sample: string, variant: AlignmentVariant) {
  const result = new Map<string, number>();
  if (variant === "uncollapsed") {
    for (const record of bundle.records) if (record.sample === sample && record.alignedNt) {
      result.set(record.id, 1);
      result.set(uncollapsedSequenceName(record), 1);
    }
  } else {
    for (const group of bundle.collapseGroups?.[sample] ?? []) {
      result.set(group.representativeId, group.familyCount);
      if (group.functionalPass) result.set(functionalSequenceName(group), group.familyCount);
    }
  }
  if (variant === "functional") {
    const reference = bundle.referenceAlignments?.[`${sample}/functional-nucleotide`];
    if (reference) {
      try {
        const row = inspectAlignment(reference, 1).records[0];
        if (row) result.set(row.name, 0);
      } catch { /* malformed result files are rejected before export */ }
    }
  }
  return result;
}

function modalTip(records: readonly { name: string; sequence: string }[], safeNames: readonly string[], counts: ReadonlyMap<string, number>) {
  const groups = new Map<string, { weight: number; index: number }>();
  records.forEach((record, index) => {
    const sequence = record.sequence.toUpperCase(), weight = counts.get(record.name) ?? 1, current = groups.get(sequence);
    if (current) { if (current.weight === 0 && weight > 0) current.index = index; current.weight += weight; }
    else groups.set(sequence, { weight, index });
  });
  let best: { sequence: string; weight: number; index: number } | undefined;
  for (const [sequence, group] of groups) if (!best || group.weight > best.weight) best = { sequence, ...group };
  return best ? { ...best, treeName: safeNames[best.index] } : undefined;
}

function empty(title: string, detail: string) {
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${WIDTH} 260" role="img"><title>${xml(title)}</title><rect width="100%" height="100%" fill="#ffffff"/><text x="40" y="80" font-family="Arial,sans-serif" font-size="24" font-weight="700">${xml(title)}</text><text x="40" y="125" font-family="Arial,sans-serif" font-size="18">${xml(detail)}</text></svg>`;
}

/**
 * Publication-style static phylogram + modal highlighter. The tree and matrix
 * each receive one half of the fixed canvas, independently of alignment width.
 * Every alignment row is reordered to its corresponding tree leaf.
 */
export function staticTreeHighlighterSvg(bundle: ResultBundle, sample: string, variant: AlignmentVariant): string {
  const active = effectiveAlignment(bundle, sample, variant), label = variant === "collapsed" ? "collapsed" : variant === "uncollapsed" ? "uncollapsed" : "functional";
  if (!active.fasta) return empty(`${sample} ${label} tree + highlighter`, `No ${label} alignment is available.`);
  const newick = active.edit?.treeNewick ?? bundle.trees[active.key];
  if (!newick) return empty(`${sample} ${label} tree + highlighter`, `The ${label} phylogeny has not been inferred.`);

  try {
    const records = inspectAlignment(active.fasta, 1).records, safeNames = treeTipNames(records.map((record) => record.name));
    let root = parseNewick(newick);
    const expected = new Set(safeNames), actual = leaves(root);
    if (actual.length !== expected.size || actual.some((name) => !expected.has(name)))
      return empty(`${sample} ${label} tree + highlighter`, "The active alignment and stored tree contain different sequence sets. Recalculate the tree first.");

    const counts = familyCounts(bundle, sample, variant), modal = modalTip(records, safeNames, counts);
    if (modal) root = rootOnOutgroup(root, modal.treeName);
    root = ladderizeTree(root, "small-first");
    const byTreeName = new Map(safeNames.map((name, index) => [name, records[index]]));
    const layout = layoutTree(root, TREE_BRANCH_WIDTH, ROW_HEIGHT, "phylogram", 0, TOP);
    const leafNodes = layout.nodes.filter((node) => !node.children.length).sort((left, right) => left.y - right.y);
    const bottom = TOP + Math.max(0, leafNodes.length - 1) * ROW_HEIGHT, height = Math.max(270, bottom + 105);
    const title = `${sample} ${label} phylogram and modal-sequence highlighter`;
    const parts: string[] = [
      `<svg xmlns="http://www.w3.org/2000/svg" width="${WIDTH}" height="${height}" viewBox="0 0 ${WIDTH} ${height}" role="img" aria-label="${xml(title)}">`,
      `<title>${xml(title)}</title><desc>The tree and modal-sequence highlighter each occupy half the plot width. Small clades are shown first. Tip-circle area is exactly ${BUBBLE_AREA_PER_FAMILY} square pixels per represented UMI family; variants above ten percent of sampled families are red.</desc>`,
      `<rect width="${WIDTH}" height="${height}" fill="#ffffff"/>`,
    ];

    // Alignment coordinate ruler: zero-based alignment columns, matching the
    // static report convention in the supplied reference figure.
    const columns = records[0]?.sequence.length ?? 0, coordinateStep = niceStep(columns);
    parts.push(`<line x1="${MATRIX_LEFT}" x2="${MATRIX_LEFT + MATRIX_WIDTH}" y1="44" y2="44" stroke="#111" stroke-width="3"/>`);
    for (let value = 0; value <= columns; value += coordinateStep) {
      const x = MATRIX_LEFT + value / Math.max(1, columns) * MATRIX_WIDTH;
      parts.push(`<line x1="${x}" x2="${x}" y1="44" y2="57" stroke="#111" stroke-width="3"/><text x="${x}" y="33" text-anchor="middle" font-family="Arial,sans-serif" font-size="18" font-weight="700">${value}</text>`);
      const minor = value + coordinateStep / 2;
      if (minor < columns) { const minorX = MATRIX_LEFT + minor / columns * MATRIX_WIDTH; parts.push(`<line x1="${minorX}" x2="${minorX}" y1="44" y2="55" stroke="#111" stroke-width="2"/>`); }
    }
    if (columns % coordinateStep) parts.push(`<line x1="${MATRIX_LEFT + MATRIX_WIDTH}" x2="${MATRIX_LEFT + MATRIX_WIDTH}" y1="44" y2="57" stroke="#111" stroke-width="3"/>`);

    // One horizontal guide per leaf makes the exact tree/alignment row mapping
    // visually explicit even where a row matches the modal sequence entirely.
    for (const node of leafNodes) parts.push(`<line x1="${MATRIX_LEFT}" x2="${MATRIX_LEFT + MATRIX_WIDTH}" y1="${node.y}" y2="${node.y}" stroke="#bdbdbd" stroke-width="0.8"/>`);
    for (const edge of layout.edges)
      parts.push(`<path d="M${edge.parent.x},${edge.parent.y} V${edge.child.y} H${edge.child.x}" fill="none" stroke="#111" stroke-width="1.25"/>`);

    const modalSequence = modal?.sequence ?? records[0]?.sequence.toUpperCase() ?? "";
    const representedFamilyTotal = records.reduce((sum, record) => sum + Math.max(0, counts.get(record.name) ?? 1), 0);
    for (const node of leafNodes) {
      const record = byTreeName.get(node.name);
      if (!record) continue;
      const count = counts.get(record.name) ?? 1, isReference = count === 0;
      if (isReference) {
        parts.push(`<rect x="${node.x - 3.5}" y="${node.y - 3.5}" width="7" height="7" fill="#ffffff" stroke="#d49a19" stroke-width="1.5"><title>Functional reference</title></rect>`);
      } else {
        const radius = Math.sqrt(count * BUBBLE_AREA_PER_FAMILY / Math.PI), abundant = count / Math.max(1, representedFamilyTotal) > .1;
        const fill = abundant ? "#FF0000" : "#000000";
        parts.push(`<circle cx="${node.x}" cy="${node.y}" r="${radius}" fill="${fill}" fill-opacity="0.6"><title>${xml(`${record.name}; ${count} UMI ${count === 1 ? "family" : "families"}`)}</title></circle>`);
      }
      const abundant = !isReference && count / Math.max(1, representedFamilyTotal) > .1;
      parts.push(`<text x="${node.x + (isReference ? 7 : Math.sqrt(Math.max(0, count) * BUBBLE_AREA_PER_FAMILY / Math.PI) + 4)}" y="${node.y + 3}" font-family="Arial,sans-serif" font-size="8" font-weight="600" fill="${abundant ? "#FF0000" : "#111"}">${xml(record.name)}</text>`);

      const sequence = record.sequence.toUpperCase();
      let start = 0;
      const fillAt = (column: number) => {
        const base = sequence[column] ?? "-";
        if (base === (modalSequence[column] ?? "-").toUpperCase()) return "";
        return COLORS[base] ?? "#737373";
      };
      while (start < columns) {
        const fill = fillAt(start); let end = start + 1;
        while (end < columns && fillAt(end) === fill) end++;
        if (fill) {
          const x = MATRIX_LEFT + start / columns * MATRIX_WIDTH, width = (end - start) / columns * MATRIX_WIDTH;
          parts.push(`<rect x="${x}" y="${node.y - ROW_HEIGHT * .42}" width="${Math.max(width, .18)}" height="${ROW_HEIGHT * .84}" fill="${fill}"/>`);
        }
        start = end;
      }
    }

    const scale = niceStep(layout.maximumDistance), scalePixels = scale / Math.max(layout.maximumDistance, scale) * (TREE_BRANCH_WIDTH - 48), footerY = bottom + 55;
    if (layout.maximumDistance > 0) parts.push(`<line x1="320" x2="${320 + scalePixels}" y1="${footerY}" y2="${footerY}" stroke="#111" stroke-width="4"/><text x="${320 + scalePixels / 2}" y="${footerY + 24}" text-anchor="middle" font-family="Arial,sans-serif" font-size="14">${scale.toPrecision(2)} substitutions/site</text>`);
    parts.push(`<text x="20" y="${height - 24}" font-family="Arial,sans-serif" font-size="11">Tip-circle area = ${BUBBLE_AREA_PER_FAMILY} px² × represented UMI-family count (no floor or cap)</text>`);
    parts.push(`<text x="20" y="${height - 9}" font-family="Arial,sans-serif" font-size="10" fill="#FF0000">Red = variant represents &gt;10% of sampled UMI families · tip order = small clades first</text>`);
    const legend = [["A", COLORS.A], ["G", COLORS.G], ["T", COLORS.T], ["C", COLORS.C], ["−", COLORS["-"]]] as const;
    legend.forEach(([base, color], index) => {
      const x = MATRIX_LEFT + 25 + index * 140;
      parts.push(`<text x="${x}" y="${footerY + 7}" font-family="serif" font-size="20">${base}</text><rect x="${x + 28}" y="${footerY - 17}" width="31" height="31" fill="${color}"/>`);
    });
    parts.push(`<line x1="${TREE_PANE_WIDTH}" x2="${TREE_PANE_WIDTH}" y1="20" y2="${bottom + 18}" stroke="#eeeeee"/>`, "</svg>");
    return parts.join("");
  } catch (cause) {
    return empty(`${sample} ${label} tree + highlighter`, `The static figure could not be generated: ${cause instanceof Error ? cause.message : String(cause)}`);
  }
}
