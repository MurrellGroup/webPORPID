import { useEffect, useMemo, useRef, useState } from "react";
import { modalSequence, referenceCoordinateLabels, referenceDisplayColumns, type CoordinateUnits } from "../alignment-regions";
import { translateAlignedNucleotides, type AlignmentFrameOffset } from "../alignment-utils";
import { parseFasta } from "../config";
import { aminoAcidBranchMutations, cladeSignature, mapParsimonyMutations, type BranchMutation } from "../phylo-mutations";
import { ladderizeTree, layoutTree, midpointRoot, parseNewick, rootOnOutgroup, serializeNewick, type LadderizeDirection, type TreeLayoutMode, type TreeNode } from "../tree";
import { treeTipNames } from "../tree-names";

export interface LeafMetadata {
  familyCount?: number;
  minimumAgreement?: number;
  /** Optional categorical color, used by contamination-reference phylogenies. */
  color?: string;
  category?: string;
}

export interface TipLegendEntry { label: string; color: string }

export interface ModalRootTip {
  treeName: string;
  sequenceName: string;
  modalSequence: string;
  representedFamilies: number;
  matchingTips: number;
}

const NT_COLORS: Record<string, string> = { A: "#78c679", C: "#6baed6", G: "#f2c14e", T: "#ef7668", U: "#ef7668", N: "#c7cbc7", "-": "#ffffff" };
const AA_GROUPS: Array<[string, string]> = [["AILMFWV", "#79b8a4"], ["KRH", "#6baed6"], ["DE", "#ef7668"], ["STNQ", "#c9df75"], ["CGP", "#e5b85c"], ["Y", "#a98dcc"], ["X*-", "#c7cbc7"]];
const RULER_HEIGHT = 30;
const DEFAULT_BUBBLE_AREA_PER_FAMILY = 25;
const aaColor = (base: string) => AA_GROUPS.find(([letters]) => letters.includes(base.toUpperCase()))?.[1] ?? "#ddd8cc";
const shortName = (name: string, maximum = 28) => name.length <= maximum ? name : `${name.slice(0, maximum - 1)}…`;

function downloadText(name: string, text: string, mime: string) {
  const url = URL.createObjectURL(new Blob([text], { type: mime })), anchor = document.createElement("a");
  anchor.href = url; anchor.download = name; anchor.click(); setTimeout(() => URL.revokeObjectURL(url), 0);
}

function saveSvg(svg: SVGSVGElement | null, name: string) {
  if (!svg) return;
  const clone = svg.cloneNode(true) as SVGSVGElement; clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  downloadText(name, new XMLSerializer().serializeToString(clone), "image/svg+xml;charset=utf-8");
}

/** Circle area is exactly familyCount × areaPerFamily; no visual cap or floor. */
export function abundanceBubbleRadius(familyCount: number, areaPerFamily = DEFAULT_BUBBLE_AREA_PER_FAMILY): number {
  return Math.sqrt(familyCount * areaPerFamily / Math.PI);
}

/** Choose a real tip carrying the most frequent observed aligned sequence. */
export function chooseModalRootTip(
  records: readonly { name: string; sequence: string }[],
  treeNames: readonly string[],
  metadata: Readonly<Record<string, LeafMetadata>> = {},
): ModalRootTip | undefined {
  if (!records.length || records.length !== treeNames.length) return undefined;
  const groups = new Map<string, { firstIndex: number; representedFamilies: number; matchingTips: number }>();
  records.forEach((record, index) => {
    const sequence = record.sequence.toUpperCase();
    const rawWeight = metadata[record.name]?.familyCount ?? 1;
    const weight = Number.isFinite(rawWeight) && rawWeight > 0 ? rawWeight : 1;
    const group = groups.get(sequence);
    if (group) { group.representedFamilies += weight; group.matchingTips += 1; }
    else groups.set(sequence, { firstIndex: index, representedFamilies: weight, matchingTips: 1 });
  });
  // In a tie, retain the earliest FASTA row. This is deterministic and does
  // not depend on rerooted traversal or ladderization order.
  const [modal, best] = [...groups].reduce((winner, candidate) =>
    candidate[1].representedFamilies > winner[1].representedFamilies ? candidate : winner);
  const bestIndex = best.firstIndex;
  return {
    treeName: treeNames[bestIndex], sequenceName: records[bestIndex].name, modalSequence: modal,
    representedFamilies: best.representedFamilies, matchingTips: best.matchingTips,
  };
}

const xml = (value: unknown) => String(value).replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll(">", "&gt;").replaceAll('"', "&quot;");

function cellTextColor(background: string): string {
  if (!background.startsWith("#") || background.length !== 7) return "#132321";
  const red = Number.parseInt(background.slice(1, 3), 16), green = Number.parseInt(background.slice(3, 5), 16), blue = Number.parseInt(background.slice(5, 7), 16);
  return red * .299 + green * .587 + blue * .114 < 135 ? "#fff" : "#132321";
}

function leaves(root: TreeNode): string[] {
  const output: string[] = [];
  const visit = (node: TreeNode) => node.children.length ? node.children.forEach(visit) : output.push(node.name);
  visit(root); return output;
}

function useViewportSize(ref: React.RefObject<HTMLDivElement | null>) {
  const [size, setSize] = useState({ width: 1, height: 1 });
  useEffect(() => {
    const node = ref.current; if (!node) return;
    const update = () => setSize({ width: Math.max(1, node.clientWidth), height: Math.max(1, node.clientHeight) }); update();
    const observer = new ResizeObserver(update); observer.observe(node); return () => observer.disconnect();
  }, [ref]); return size;
}

function canvasResolution(canvas: HTMLCanvasElement, width: number, height: number) {
  const ratio = window.devicePixelRatio || 1;
  canvas.width = Math.floor(width * ratio); canvas.height = Math.floor(height * ratio);
  canvas.style.width = `${width}px`; canvas.style.height = `${height}px`;
  const context = canvas.getContext("2d")!; context.setTransform(ratio, 0, 0, ratio, 0, 0); return context;
}

function AlignmentCanvas({ sequences, columns, labels, modal, alphabet, highlighter, rowHeight, cellWidth, selected, onSelect, scrollRef, peerRef }: {
  sequences: string[]; columns: number[]; labels: string[]; modal: string; alphabet: "nt" | "aa"; highlighter: boolean;
  rowHeight: number; cellWidth: number; selected: number; onSelect(index: number): void;
  scrollRef: React.RefObject<HTMLDivElement | null>; peerRef: React.RefObject<HTMLDivElement | null>;
}) {
  const canvas = useRef<HTMLCanvasElement>(null), size = useViewportSize(scrollRef);
  useEffect(() => {
    const container = scrollRef.current, target = canvas.current; if (!container || !target) return;
    const draw = () => {
      const context = canvasResolution(target, size.width, size.height), left = container.scrollLeft, top = container.scrollTop;
      const firstRow = Math.max(0, Math.floor((top - RULER_HEIGHT) / rowHeight)), lastRow = Math.min(sequences.length, Math.ceil((top + size.height - RULER_HEIGHT) / rowHeight) + 1);
      const firstDisplayColumn = Math.max(0, Math.floor(left / cellWidth)), lastDisplayColumn = Math.min(columns.length, Math.ceil((left + size.width) / cellWidth) + 1);
      context.clearRect(0, 0, size.width, size.height); context.textAlign = "center"; context.textBaseline = "middle";
      context.font = `700 ${Math.min(11, rowHeight - 4, cellWidth - 2)}px ui-monospace,monospace`;
      for (let row = firstRow; row < lastRow; row += 1) for (let display = firstDisplayColumn; display < lastDisplayColumn; display += 1) {
        const column = columns[display], base = sequences[row]?.[column] ?? "-", x = display * cellWidth - left, y = RULER_HEIGHT + row * rowHeight - top;
        const same = base.toUpperCase() === (modal[column] ?? "-").toUpperCase();
        context.fillStyle = highlighter && same ? (base === "-" ? "#ffffff" : "#eceeea") : alphabet === "nt" ? NT_COLORS[base.toUpperCase()] ?? "#c7cbc7" : aaColor(base);
        context.fillRect(x, y, cellWidth + .4, rowHeight + .4);
        if (cellWidth >= 9) { context.fillStyle = highlighter && same ? "#8a938f" : "#132321"; context.fillText(base, x + cellWidth / 2, y + rowHeight / 2); }
      }
      if (selected >= firstRow && selected < lastRow) {
        context.strokeStyle = "#132321"; context.lineWidth = 2; context.strokeRect(0, RULER_HEIGHT + selected * rowHeight - top, size.width, rowHeight);
      }
      // The canvas is viewport-sticky. Draw the reference ruler last at a
      // fixed viewport position so sequence rows cannot paint over it while
      // the user scrolls vertically.
      context.fillStyle = "#f1f2ed"; context.fillRect(0, 0, size.width, RULER_HEIGHT);
      context.font = "8px ui-monospace,monospace"; context.fillStyle = "#596862";
      const tickEvery = Math.max(1, Math.ceil(44 / cellWidth));
      for (let display = firstDisplayColumn; display < lastDisplayColumn; display += 1) {
        const original = columns[display], x = display * cellWidth - left;
        if (display === 0 || display === columns.length - 1 || display % tickEvery === 0 || (display && original !== columns[display - 1] + 1)) {
          context.fillText(labels[original] || String(original + 1), x + cellWidth / 2, 10);
          context.strokeStyle = "#aeb8b2"; context.beginPath(); context.moveTo(x + cellWidth / 2, 18); context.lineTo(x + cellWidth / 2, RULER_HEIGHT); context.stroke();
        }
      }
    };
    draw();
    const sync = () => {
      if (peerRef.current && Math.abs(peerRef.current.scrollTop - container.scrollTop) > 1) peerRef.current.scrollTop = container.scrollTop;
      draw();
    };
    container.addEventListener("scroll", sync, { passive: true }); return () => container.removeEventListener("scroll", sync);
  }, [sequences, columns, labels, modal, alphabet, highlighter, rowHeight, cellWidth, selected, size, scrollRef, peerRef]);
  return <div ref={scrollRef} className="alignment-viewport" onPointerDown={(event) => {
    const index = Math.floor((scrollRef.current!.scrollTop + event.nativeEvent.offsetY - RULER_HEIGHT) / rowHeight);
    if (index >= 0 && index < sequences.length) onSelect(index);
  }}><div style={{ width: Math.max(size.width, columns.length * cellWidth), height: Math.max(size.height, RULER_HEIGHT + sequences.length * rowHeight), position: "relative" }}><canvas ref={canvas} /></div></div>;
}

function niceScale(maximum: number): number {
  if (!(maximum > 0)) return 0;
  const target = maximum / 4, power = 10 ** Math.floor(Math.log10(target)), scaled = target / power;
  return (scaled >= 5 ? 5 : scaled >= 2 ? 2 : 1) * power;
}

function agreementColor(value: number): string {
  const amount = Math.max(0, Math.min(1, (value - .5) / .5));
  return `rgb(${Math.round(222 - 190 * amount)},${Math.round(105 + 55 * amount)},${Math.round(82 + 45 * amount)})`;
}

function TreeSvg({ root, width, rowHeight, selected, onSelect, scrollRef, peerRef, svgRef, layoutMode, showNames, leafMetadata, leafLabels, bubbleAreaPerFamily, colorByAgreement, mutations, mutationLimit }: {
  root: TreeNode; width: number; rowHeight: number; selected: number; onSelect(index: number): void;
  scrollRef: React.RefObject<HTMLDivElement | null>; peerRef: React.RefObject<HTMLDivElement | null>; svgRef: React.RefObject<SVGSVGElement | null>;
  layoutMode: TreeLayoutMode; showNames: boolean; leafMetadata: Readonly<Record<string, LeafMetadata>>; leafLabels: Readonly<Record<string, string>>; bubbleAreaPerFamily: number; colorByAgreement: boolean;
  mutations: ReadonlyMap<string, BranchMutation[]>; mutationLimit: number;
}) {
  const labelWidth = showNames ? Math.min(220, Math.max(120, width * .36)) : 14, branchWidth = Math.max(80, width - labelWidth);
  const layout = useMemo(() => layoutTree(root, branchWidth, rowHeight, layoutMode, 18, RULER_HEIGHT + rowHeight / 2), [root, branchWidth, rowHeight, layoutMode]);
  const leafNodes = useMemo(() => layout.nodes.filter((node) => !node.children.length).sort((left, right) => left.y - right.y), [layout]);
  const scale = niceScale(layout.maximumDistance), horizontalPadding = Math.min(24, branchWidth * .1);
  const scalePixels = scale && layout.maximumDistance ? scale / layout.maximumDistance * Math.max(0, branchWidth - 2 * horizontalPadding) : 0;
  useEffect(() => {
    const container = scrollRef.current; if (!container) return;
    const sync = () => { if (peerRef.current && Math.abs(peerRef.current.scrollTop - container.scrollTop) > 1) peerRef.current.scrollTop = container.scrollTop; };
    container.addEventListener("scroll", sync, { passive: true }); return () => container.removeEventListener("scroll", sync);
  }, [scrollRef, peerRef]);
  const height = Math.max(layout.height, RULER_HEIGHT + leafNodes.length * rowHeight);
  const abundanceScaleApplies = Object.values(leafMetadata).some((metadata) => (metadata.familyCount ?? 1) > 1);
  return <div ref={scrollRef} className="tree-viewport" onPointerDown={(event) => {
    const index = Math.floor((scrollRef.current!.scrollTop + event.nativeEvent.offsetY - RULER_HEIGHT) / rowHeight);
    if (index >= 0 && index < leafNodes.length) onSelect(index);
  }}><svg ref={svgRef} width={width} height={height} viewBox={`0 0 ${width} ${height}`} role="img" aria-label="Phylogram coordinated with aligned sequences">
    <desc>Tip circle area is family count multiplied by {bubbleAreaPerFamily} square pixels, without a radius floor, cap, or dataset normalization.</desc>
    <rect width={width} height={height} fill="#fbfaf5" />
    {abundanceScaleApplies && <text x={width - 6} y="12" textAnchor="end" fontFamily="ui-monospace,monospace" fontSize="7" fill="#466047">area = {bubbleAreaPerFamily} px² × UMI families</text>}
    {selected >= 0 && <rect x="0" y={RULER_HEIGHT + selected * rowHeight} width={width} height={rowHeight} fill="rgba(201,239,115,.32)" />}
    {scalePixels >= 8 && <g><line x1={horizontalPadding} x2={horizontalPadding + scalePixels} y1="11" y2="11" stroke="#273d37" strokeWidth="1.2" /><line x1={horizontalPadding} x2={horizontalPadding} y1="7" y2="15" stroke="#273d37" /><line x1={horizontalPadding + scalePixels} x2={horizontalPadding + scalePixels} y1="7" y2="15" stroke="#273d37" /><text x={horizontalPadding + scalePixels / 2} y="25" textAnchor="middle" fontFamily="Inter,Arial,sans-serif" fontSize="8" fill="#4e5d58">{scale.toPrecision(2)} substitutions/site</text></g>}
    <line x1={branchWidth} x2={branchWidth} y1={RULER_HEIGHT} y2={height} stroke="#d2d8d3" />
    {layout.edges.map((edge, index) => {
      const branch = mutations.get(cladeSignature(edge.child)) ?? [], shown = branch.slice(0, mutationLimit);
      const label = shown.map((mutation) => `${mutation.column + 1}${mutation.from}→${mutation.to}`).join(" · ") + (branch.length > mutationLimit ? ` · +${branch.length - mutationLimit}` : "");
      const x = edge.parent.x + (edge.child.x - edge.parent.x) * .5;
      return <g key={`edge-${index}`}><path d={`M${edge.parent.x},${edge.parent.y} V${edge.child.y} H${edge.child.x}`} stroke="#3f5650" strokeWidth=".95" fill="none" />{label && <text x={x} y={edge.child.y - 5} textAnchor="middle" fontFamily="ui-monospace,monospace" fontSize="7" fontWeight="650" fill="#71374a" stroke="#fbfaf5" strokeWidth="2.5" paintOrder="stroke fill">{label}</text>}</g>;
    })}
    {layout.nodes.map((node, index) => {
      const leaf = !node.children.length, metadata = leafMetadata[node.name] ?? {}, count = metadata.familyCount ?? 1, label = leafLabels[node.name] ?? node.name;
      const radius = leaf ? abundanceBubbleRadius(count, bubbleAreaPerFamily) : 1.6;
      const fill = metadata.color ?? (colorByAgreement && metadata.minimumAgreement != null ? agreementColor(metadata.minimumAgreement) : "#08796f");
      return <g key={`node-${index}`}><circle cx={node.x} cy={node.y} r={radius} fill={leaf ? fill : "#fbfaf5"} fillOpacity=".9" stroke="#31534b" strokeWidth=".8"><title>{leaf ? `${label} · ${count} UMI ${count === 1 ? "family" : "families"}${metadata.minimumAgreement == null ? "" : ` · family minimum agreement ${metadata.minimumAgreement.toFixed(3)}`}${metadata.category ? ` · ${metadata.category}` : ""}` : `${cladeSignature(node).split("\u0000").length} descendants`}</title></circle>{leaf && showNames && <><line x1={node.x + radius + 2} x2={branchWidth + 4} y1={node.y} y2={node.y} stroke="#b3bdb8" strokeWidth=".65" strokeDasharray="2 3" /><text x={branchWidth + 8} y={node.y + 3.3} fontFamily="ui-monospace,monospace" fontSize="9" fill="#263630">{shortName(label)}<title>{label}</title></text></>}</g>;
    })}
  </svg></div>;
}

interface TreeAlignmentSvgOptions {
  root: TreeNode;
  treeWidth: number;
  rowHeight: number;
  cellWidth: number;
  layoutMode: TreeLayoutMode;
  showNames: boolean;
  leafMetadata: Readonly<Record<string, LeafMetadata>>;
  leafLabels: Readonly<Record<string, string>>;
  sequencesByTip: Readonly<Record<string, string>>;
  columns: readonly number[];
  labels: readonly string[];
  modal: string;
  alphabet: "nt" | "aa";
  highlighter: boolean;
  bubbleAreaPerFamily: number;
  showAbundanceScale: boolean;
  colorByAgreement: boolean;
  mutations: ReadonlyMap<string, BranchMutation[]>;
  mutationLimit: number;
  tipLegend: readonly TipLegendEntry[];
}

/** Build the same coordinated tree-and-residue SVG used on screen, without rasterizing the canvas alignment. */
export function treeAlignmentSvg(options: TreeAlignmentSvgOptions): string {
  const { root, treeWidth, rowHeight, cellWidth, layoutMode, showNames, leafMetadata, leafLabels, sequencesByTip,
    columns, labels, modal, alphabet, highlighter, bubbleAreaPerFamily, showAbundanceScale, colorByAgreement, mutations, mutationLimit, tipLegend } = options;
  const labelWidth = showNames ? Math.min(220, Math.max(120, treeWidth * .36)) : 14;
  const branchWidth = Math.max(80, treeWidth - labelWidth), matrixX = treeWidth + 8, matrixWidth = Math.max(cellWidth, columns.length * cellWidth);
  const width = matrixX + matrixWidth + 12, categoryLegendHeight = tipLegend.length ? 8 + Math.ceil(tipLegend.length / 4) * 18 : 0;
  const abundanceLegendHeight = showAbundanceScale ? 25 : 0, legendHeight = abundanceLegendHeight + categoryLegendHeight;
  const layout = layoutTree(root, branchWidth, rowHeight, layoutMode, 18, RULER_HEIGHT + rowHeight / 2);
  const height = Math.max(layout.height, RULER_HEIGHT + layout.leaves * rowHeight) + legendHeight;
  const scale = niceScale(layout.maximumDistance), horizontalPadding = Math.min(24, branchWidth * .1);
  const scalePixels = scale && layout.maximumDistance ? scale / layout.maximumDistance * Math.max(0, branchWidth - 2 * horizontalPadding) : 0;
  const fragments: string[] = [`<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" role="img" aria-label="Tree coordinated with aligned leaf sequences">`, `<rect width="${width}" height="${height}" fill="#fbfaf5"/>`];
  if (scalePixels >= 8) fragments.push(`<g><line x1="${horizontalPadding}" x2="${horizontalPadding + scalePixels}" y1="11" y2="11" stroke="#273d37" stroke-width="1.2"/><line x1="${horizontalPadding}" x2="${horizontalPadding}" y1="7" y2="15" stroke="#273d37"/><line x1="${horizontalPadding + scalePixels}" x2="${horizontalPadding + scalePixels}" y1="7" y2="15" stroke="#273d37"/><text x="${horizontalPadding + scalePixels / 2}" y="25" text-anchor="middle" font-family="Inter,Arial,sans-serif" font-size="8" fill="#4e5d58">${xml(scale.toPrecision(2))} substitutions/site</text></g>`);
  const tickEvery = Math.max(1, Math.ceil(44 / cellWidth));
  columns.forEach((column, display) => {
    if (display === 0 || display === columns.length - 1 || display % tickEvery === 0 || (display && column !== columns[display - 1] + 1)) {
      const x = matrixX + display * cellWidth + cellWidth / 2;
      fragments.push(`<text transform="translate(${x},25) rotate(-55)" text-anchor="end" font-family="ui-monospace,monospace" font-size="7" fill="#596862">${xml(labels[column] || column + 1)}</text>`);
    }
  });
  fragments.push(`<line x1="${branchWidth}" x2="${branchWidth}" y1="${RULER_HEIGHT}" y2="${height - legendHeight}" stroke="#d2d8d3"/>`);
  layout.edges.forEach((edge) => {
    const branch = mutations.get(cladeSignature(edge.child)) ?? [], shown = branch.slice(0, mutationLimit);
    const label = shown.map((mutation) => `${mutation.column + 1}${mutation.from}→${mutation.to}`).join(" · ") + (branch.length > mutationLimit ? ` · +${branch.length - mutationLimit}` : "");
    const x = edge.parent.x + (edge.child.x - edge.parent.x) * .5;
    fragments.push(`<path d="M${edge.parent.x},${edge.parent.y} V${edge.child.y} H${edge.child.x}" stroke="#3f5650" stroke-width=".95" fill="none"/>`);
    if (label) fragments.push(`<text x="${x}" y="${edge.child.y - 5}" text-anchor="middle" font-family="ui-monospace,monospace" font-size="7" font-weight="650" fill="#71374a" stroke="#fbfaf5" stroke-width="2.5" paint-order="stroke fill">${xml(label)}</text>`);
  });
  layout.nodes.forEach((node) => {
    const leaf = !node.children.length, metadata = leafMetadata[node.name] ?? {}, count = metadata.familyCount ?? 1;
    const radius = leaf ? abundanceBubbleRadius(count, bubbleAreaPerFamily) : 1.6, label = leafLabels[node.name] ?? node.name;
    const fill = metadata.color ?? (colorByAgreement && metadata.minimumAgreement != null ? agreementColor(metadata.minimumAgreement) : "#08796f");
    if (leaf && showNames) fragments.push(`<line x1="${node.x + radius + 2}" x2="${branchWidth + 4}" y1="${node.y}" y2="${node.y}" stroke="#b3bdb8" stroke-width=".65" stroke-dasharray="2 3"/><text x="${branchWidth + 8}" y="${node.y + 3.3}" font-family="ui-monospace,monospace" font-size="9" fill="#263630">${xml(shortName(label))}<title>${xml(label)}</title></text>`);
    fragments.push(`<circle cx="${node.x}" cy="${node.y}" r="${radius}" fill="${leaf ? xml(fill) : "#fbfaf5"}" fill-opacity=".9" stroke="#31534b" stroke-width=".8"><title>${xml(leaf ? `${label} · ${count} UMI ${count === 1 ? "family" : "families"}${metadata.category ? ` · ${metadata.category}` : ""}` : `${cladeSignature(node).split("\u0000").length} descendants`)}</title></circle>`);
    if (!leaf) return;
    const sequence = sequencesByTip[node.name] ?? "";
    columns.forEach((column, display) => {
      const value = sequence[column] ?? "-", same = value.toUpperCase() === (modal[column] ?? "-").toUpperCase();
      const background = highlighter && same ? (value === "-" ? "#ffffff" : "#eceeea") : alphabet === "nt" ? NT_COLORS[value.toUpperCase()] ?? "#c7cbc7" : aaColor(value);
      const x = matrixX + display * cellWidth, y = node.y - rowHeight / 2;
      fragments.push(`<rect x="${x}" y="${y}" width="${cellWidth}" height="${rowHeight}" fill="${background}"/>`);
      if (cellWidth >= 9) fragments.push(`<text x="${x + cellWidth / 2}" y="${node.y + Math.min(4, rowHeight * .16)}" text-anchor="middle" font-family="ui-monospace,monospace" font-size="${Math.min(10.5, cellWidth - 2)}" font-weight="700" fill="${highlighter && same ? "#8a938f" : cellTextColor(background)}">${xml(value)}</text>`);
    });
  });
  const legendTop = height - legendHeight;
  if (showAbundanceScale) fragments.push(`<line x1="12" x2="${width - 12}" y1="${legendTop + 3}" y2="${legendTop + 3}" stroke="#c9d0ca"/><text x="14" y="${legendTop + 18}" font-family="Inter,Arial,sans-serif" font-size="8" font-weight="700" fill="#466047">Exact tip area = ${bubbleAreaPerFamily} px² × UMI-family count; no radius floor, cap, or dataset normalization</text>`);
  if (tipLegend.length) {
    const top = legendTop + abundanceLegendHeight, cell = (width - 24) / Math.min(4, tipLegend.length);
    tipLegend.forEach((entry, index) => {
      const column = index % 4, row = Math.floor(index / 4), x = 18 + column * cell, y = top + 14 + row * 18;
      fragments.push(`<circle cx="${x}" cy="${y - 3}" r="4" fill="${xml(entry.color)}" stroke="#244d45" stroke-width=".6"/><text x="${x + 9}" y="${y}" font-family="Inter,Arial,sans-serif" font-size="8" fill="#34453f">${xml(entry.label)}</text>`);
    });
  }
  fragments.push("</svg>"); return fragments.join("");
}

export function AlignmentTreeViewer({ fasta, newick, alphabet = "nt", onAlphabetChange, name = "webporpid-tree", frameOffset = 0, referenceSequence, leafMetadata = {}, collapsed = false, treeStale = false, tipLegend = [], variantLabel }: {
  fasta: string; newick?: string; alphabet?: "nt" | "aa"; name?: string; frameOffset?: AlignmentFrameOffset; referenceSequence?: string;
  onAlphabetChange?(alphabet: "nt" | "aa"): void; leafMetadata?: Readonly<Record<string, LeafMetadata>>; collapsed?: boolean; treeStale?: boolean;
  tipLegend?: readonly TipLegendEntry[]; variantLabel?: string;
}) {
  const nucleotideRecords = useMemo(() => parseFasta(fasta), [fasta]);
  const displayedRecords = useMemo(() => nucleotideRecords.map((record) => ({ ...record, sequence: alphabet === "nt" ? record.sequence : translateAlignedNucleotides(record.sequence, frameOffset) })), [nucleotideRecords, alphabet, frameOffset]);
  const [cellWidth, setCellWidth] = useState(11), [rowHeight, setRowHeight] = useState(20), [treeWidth, setTreeWidth] = useState(600);
  const [selected, setSelected] = useState(0), [layoutMode, setLayoutMode] = useState<TreeLayoutMode>("phylogram"), [ladderization, setLadderization] = useState<"none" | LadderizeDirection>("none");
  const [rootMode, setRootMode] = useState<"modal" | "midpoint">("modal");
  const [highlighter, setHighlighter] = useState(false), [showNames, setShowNames] = useState(true), [showMutations, setShowMutations] = useState(false);
  const [mutationLimit, setMutationLimit] = useState(2), [bubbleAreaPerFamily, setBubbleAreaPerFamily] = useState(DEFAULT_BUBBLE_AREA_PER_FAMILY), [colorByAgreement, setColorByAgreement] = useState(false);
  const [regionDraft, setRegionDraft] = useState(""), [regionText, setRegionText] = useState(""), [coordinateUnits, setCoordinateUnits] = useState<CoordinateUnits>(alphabet);
  const treeScroll = useRef<HTMLDivElement>(null), alignmentScroll = useRef<HTMLDivElement>(null), svgRef = useRef<SVGSVGElement>(null);
  const parsed = useMemo(() => {
    if (!newick) return { tree: undefined, error: "" };
    try { return { tree: parseNewick(newick), error: "" }; }
    catch (cause) { return { tree: undefined, error: cause instanceof Error ? cause.message : String(cause) }; }
  }, [newick]);
  const safeNames = useMemo(() => treeTipNames(nucleotideRecords.map((record) => record.name)), [nucleotideRecords]);
  const treeTipMismatch = useMemo(() => {
    if (!parsed.tree) return false;
    const expected = new Set(safeNames), actual = leaves(parsed.tree);
    return actual.length !== expected.size || actual.some((name) => !expected.has(name));
  }, [parsed.tree, safeNames]);
  const modalRootTip = useMemo(() => chooseModalRootTip(nucleotideRecords, safeNames, leafMetadata), [nucleotideRecords, safeNames, leafMetadata]);
  const rooted = useMemo(() => {
    if (!parsed.tree || treeTipMismatch) return { tree: parsed.tree, error: "" };
    try {
      return { tree: rootMode === "midpoint" ? midpointRoot(parsed.tree) : modalRootTip ? rootOnOutgroup(parsed.tree, modalRootTip.treeName) : parsed.tree, error: "" };
    } catch (cause) { return { tree: parsed.tree, error: cause instanceof Error ? cause.message : String(cause) }; }
  }, [parsed.tree, treeTipMismatch, rootMode, modalRootTip]);
  const orderedTree = useMemo(() => !rooted.tree || ladderization === "none" ? rooted.tree : ladderizeTree(rooted.tree, ladderization), [rooted.tree, ladderization]);
  const displayTree = treeTipMismatch ? undefined : orderedTree;
  const orderedRecords = useMemo(() => {
    if (!displayTree) return displayedRecords;
    const byName = new Map(displayedRecords.map((record, index) => [safeNames[index], record]));
    const found = leaves(displayTree).map((leaf) => byName.get(leaf)).filter((record): record is (typeof displayedRecords)[number] => Boolean(record));
    return found.length === displayedRecords.length ? found : displayedRecords;
  }, [displayTree, displayedRecords, safeNames]);
  const leafLabels = useMemo(() => Object.fromEntries(displayedRecords.map((record, index) => [safeNames[index], record.name])), [displayedRecords, safeNames]);
  const sequencesByTip = useMemo(() => Object.fromEntries(displayedRecords.map((record, index) => [safeNames[index], record.sequence])), [displayedRecords, safeNames]);
  const normalizedLeafMetadata = useMemo(() => Object.fromEntries(nucleotideRecords.map((record, index) => {
    const metadata = leafMetadata[record.name] ?? {};
    return [safeNames[index], collapsed ? { familyCount: metadata.familyCount } : metadata];
  })), [collapsed, leafMetadata, nucleotideRecords, safeNames]);
  const ntByName = useMemo(() => new Map(nucleotideRecords.map((record, index) => [safeNames[index], record.sequence])), [nucleotideRecords, safeNames]);
  const parsimony = useMemo(() => {
    if (!displayTree) return undefined;
    try { return mapParsimonyMutations(displayTree, ntByName); } catch { return undefined; }
  }, [displayTree, ntByName]);
  const mutations = useMemo(() => {
    if (!showMutations || !displayTree || !parsimony) return new Map<string, BranchMutation[]>();
    if (alphabet === "nt") return parsimony.mutationsByClade;
    const result = new Map<string, BranchMutation[]>();
    const layout = layoutTree(displayTree, 500, rowHeight, layoutMode, 18, RULER_HEIGHT + rowHeight / 2);
    for (const edge of layout.edges) {
      const signature = cladeSignature(edge.child), parent = parsimony.sequencesByClade.get(cladeSignature(edge.parent)), child = parsimony.sequencesByClade.get(signature);
      if (parent && child) result.set(signature, aminoAcidBranchMutations(parent, child, signature, frameOffset));
    }
    return result;
  }, [showMutations, displayTree, parsimony, alphabet, frameOffset, rowHeight, layoutMode]);
  useEffect(() => { setRootMode("modal"); }, [fasta, newick]);
  useEffect(() => { setSelected(0); }, [fasta, newick, ladderization, rootMode]);
  useEffect(() => setCoordinateUnits(alphabet), [alphabet]);
  const sequences = orderedRecords.map((record) => record.sequence);
  const modalNt = modalRootTip?.modalSequence ?? modalSequence(nucleotideRecords.map((record) => record.sequence));
  const modal = alphabet === "nt" ? modalNt : translateAlignedNucleotides(modalNt, frameOffset);
  const referenceNt = referenceSequence || modalNt;
  const labels = useMemo(() => referenceCoordinateLabels(referenceNt, alphabet, frameOffset), [referenceNt, alphabet, frameOffset]);
  const columns = useMemo(() => referenceDisplayColumns(referenceNt, regionText, coordinateUnits, alphabet, frameOffset, sequences[0]?.length ?? 0), [referenceNt, regionText, coordinateUnits, alphabet, frameOffset, sequences]);
  return <section className="alignment-tree-viewer">
    {parsed.error && <div className="viewer-status">The stored tree could not be rendered ({parsed.error}); the alignment remains available.</div>}
    {rooted.error && <div className="viewer-status warning">The requested rooting could not be applied ({rooted.error}). The stored root is shown unchanged.</div>}
    {treeTipMismatch && <div className="viewer-status warning">The edited alignment no longer has the same sequence set as this tree. The stale tree is hidden until you press Recalculate tree.</div>}
    {treeStale && <div className="viewer-status warning">The alignment has changed since this tree was inferred. Branch mutations use the edited alignment, but topology and branch lengths remain stale until you press Recalculate tree.</div>}
    <header><div><span className="section-kicker">{variantLabel ?? (collapsed ? "Collapsed haplotypes" : "Uncollapsed UMI families")}</span><h3>{alphabet === "nt" ? "Nucleotide" : "Amino-acid"} tree + alignment</h3><p>{columns.length.toLocaleString()} of {(sequences[0]?.length ?? 0).toLocaleString()} columns shown</p></div>
      <div className="viewer-controls">
        {displayTree && <><div className="rooting-controls" aria-label="Tree root"><span>Root</span><button type="button" className={rootMode === "modal" ? "active" : ""} aria-pressed={rootMode === "modal"} onClick={() => setRootMode("modal")}>Modal node</button><button type="button" className={rootMode === "midpoint" ? "active" : ""} aria-pressed={rootMode === "midpoint"} onClick={() => setRootMode("midpoint")}>Midpoint root</button></div><label>Tree width <input type="range" min="260" max="1100" step="10" value={treeWidth} onChange={(event) => setTreeWidth(Number(event.target.value))} /></label><label>Tree <select value={layoutMode} onChange={(event) => setLayoutMode(event.target.value as TreeLayoutMode)}><option value="phylogram">Branch lengths</option><option value="cladogram">Topology</option></select></label><label>Tip order <select value={ladderization} onChange={(event) => setLadderization(event.target.value as typeof ladderization)}><option value="none">Newick order</option><option value="large-first">Large clades first</option><option value="small-first">Small clades first</option></select></label></>}
        <label>Columns <input type="range" min="3" max="18" value={cellWidth} onChange={(event) => setCellWidth(Number(event.target.value))} /></label><label>Rows <input type="range" min="14" max="34" value={rowHeight} onChange={(event) => setRowHeight(Number(event.target.value))} /></label>
        {collapsed && displayTree && <label>Bubble area · {bubbleAreaPerFamily} px²/family<input type="range" min="1" max="100" step="1" value={bubbleAreaPerFamily} onChange={(event) => setBubbleAreaPerFamily(Number(event.target.value))} /></label>}
        {displayTree && <><button type="button" onClick={() => saveSvg(svgRef.current, `${name}.tree.svg`)}>Tree SVG</button><button type="button" onClick={() => downloadText(`${name}.tree-alignment.svg`, treeAlignmentSvg({ root: displayTree, treeWidth, rowHeight, cellWidth, layoutMode, showNames, leafMetadata: normalizedLeafMetadata, leafLabels, sequencesByTip, columns, labels, modal, alphabet, highlighter, bubbleAreaPerFamily, showAbundanceScale: collapsed, colorByAgreement, mutations, mutationLimit, tipLegend }), "image/svg+xml;charset=utf-8")}>Tree + alignment SVG</button><button type="button" onClick={() => downloadText(`${name}.nwk`, `${serializeNewick(displayTree)};\n`, "text/plain;charset=utf-8")}>Newick</button></>}
      </div></header>
    <div className="viewer-options"><div className="alignment-mode-toggle" aria-label="Alignment alphabet"><button type="button" className={alphabet === "nt" ? "active" : ""} disabled={!onAlphabetChange} onClick={() => onAlphabetChange?.("nt")}>Nucleotide</button><button type="button" className={alphabet === "aa" ? "active" : ""} disabled={!onAlphabetChange} onClick={() => onAlphabetChange?.("aa")}>Amino acid</button></div><label><input type="checkbox" checked={highlighter} onChange={(event) => setHighlighter(event.target.checked)} /> Highlight differences from modal sequence</label><label><input type="checkbox" checked={showNames} onChange={(event) => setShowNames(event.target.checked)} /> Show sequence names</label>{displayTree && <label><input type="checkbox" checked={showMutations} onChange={(event) => setShowMutations(event.target.checked)} /> Map branch mutations</label>}{showMutations && displayTree && <label>Mutation labels <input type="range" min="1" max="12" value={mutationLimit} onChange={(event) => setMutationLimit(Number(event.target.value))} /></label>}{!collapsed && !tipLegend.length && <label>Tip color <select value={colorByAgreement ? "agreement" : "uniform"} onChange={(event) => setColorByAgreement(event.target.value === "agreement")}><option value="uniform">Uniform</option><option value="agreement">Family minimum agreement</option></select></label>}<div className="region-input"><span>Reference-region filter</span><input value={regionDraft} onChange={(event) => setRegionDraft(event.target.value)} onKeyDown={(event) => { if (event.key === "Enter") { event.preventDefault(); setRegionText(regionDraft.trim()); } }} placeholder="100;120;200-205" /><select value={coordinateUnits} onChange={(event) => setCoordinateUnits(event.target.value as CoordinateUnits)}><option value="nt">NT coordinates</option><option value="aa">AA coordinates</option></select><button type="button" onClick={() => setRegionText(regionDraft.trim())}>Apply regions</button><button type="button" disabled={!regionText && !regionDraft} onClick={() => { setRegionDraft(""); setRegionText(""); }}>Show all</button></div></div>
    {displayTree && rootMode === "modal" && modalRootTip && <div className="viewer-status root-status">Rooted on the branch immediately above the observed modal-sequence node <strong>{modalRootTip.sequenceName}</strong>, with an exactly zero-length edge to that tip. Its aligned sequence represents {modalRootTip.representedFamilies.toLocaleString()} retained UMI {modalRootTip.representedFamilies === 1 ? "family" : "families"}{modalRootTip.matchingTips > 1 ? ` across ${modalRootTip.matchingTips.toLocaleString()} identical tree tips` : ""}.</div>}
    {displayTree && rootMode === "midpoint" && <div className="viewer-status root-status">Midpoint-rooted at half the longest tip-to-tip path; topology and pairwise patristic distances are unchanged.</div>}
    {collapsed && <div className="bubble-science-note">Exact abundance scale: circle area = {bubbleAreaPerFamily} px² × retained UMI-family count. No radius floor, cap, or dataset normalization is applied.</div>}
    {tipLegend.length > 0 && <div className="viewer-tip-legend">{tipLegend.map((entry) => <span key={entry.label}><i style={{ background: entry.color }} />{entry.label}</span>)}</div>}
    {regionText && columns.length > 0 && <div className="viewer-status region-active">Reference-region filter active: {regionText} ({coordinateUnits.toUpperCase()} coordinates). Press Show all to restore every alignment column.</div>}
    {regionText && !columns.length && <div className="viewer-status warning">No displayed columns matched those reference coordinates.</div>}
    <div className={`viewer-grid ${displayTree ? "with-tree" : "alignment-only"}`} style={displayTree ? { gridTemplateColumns: `${treeWidth}px minmax(0,1fr)` } : undefined}>
      {displayTree && <TreeSvg root={displayTree} width={treeWidth} rowHeight={rowHeight} selected={selected} onSelect={setSelected} scrollRef={treeScroll} peerRef={alignmentScroll} svgRef={svgRef} layoutMode={layoutMode} showNames={showNames} leafMetadata={normalizedLeafMetadata} leafLabels={leafLabels} bubbleAreaPerFamily={bubbleAreaPerFamily} colorByAgreement={colorByAgreement} mutations={mutations} mutationLimit={mutationLimit} />}
      <AlignmentCanvas sequences={sequences} columns={columns} labels={labels} modal={modal} alphabet={alphabet} highlighter={highlighter} rowHeight={rowHeight} cellWidth={cellWidth} selected={selected} onSelect={setSelected} scrollRef={alignmentScroll} peerRef={treeScroll} />
    </div>
    <footer><strong>{orderedRecords[selected]?.name ?? "No sequence"}</strong><span>{orderedRecords.length.toLocaleString()} rows × {(sequences[0]?.length ?? 0).toLocaleString()} columns · {highlighter ? "modal highlighter" : "residue colors"} · alignment cells are canvas-virtualized</span></footer>
  </section>;
}
