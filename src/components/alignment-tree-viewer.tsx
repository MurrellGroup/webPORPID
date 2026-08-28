import { useEffect, useMemo, useRef, useState } from "react";
import { parseFasta } from "../config";
import {
  ladderizeTree, layoutTree, parseNewick, serializeNewick,
  type LadderizeDirection, type TreeLayoutMode, type TreeNode,
} from "../tree";

const NT_COLORS: Record<string, string> = { A: "#78c679", C: "#6baed6", G: "#f2c14e", T: "#ef7668", U: "#ef7668", N: "#c7cbc7", "-": "#f7f5ef" };
const AA_GROUPS: Array<[string, string]> = [["AILMFWV", "#79b8a4"], ["KRH", "#6baed6"], ["DE", "#ef7668"], ["STNQ", "#c9df75"], ["CGP", "#e5b85c"], ["Y", "#a98dcc"], ["X*-", "#c7cbc7"]];
const aaColor = (base: string) => AA_GROUPS.find(([letters]) => letters.includes(base.toUpperCase()))?.[1] ?? "#ddd8cc";
const safeTreeName = (name: string) => name.replace(/[^A-Za-z0-9_.|*+\-]/g, "_");
const shortName = (name: string, maximum = 24) => name.length <= maximum ? name : `${name.slice(0, maximum - 1)}…`;

function downloadText(name: string, text: string, mime: string) {
  const url = URL.createObjectURL(new Blob([text], { type: mime })), anchor = document.createElement("a");
  anchor.href = url; anchor.download = name; anchor.click(); setTimeout(() => URL.revokeObjectURL(url), 0);
}

function saveSvg(svg: SVGSVGElement | null, name: string) {
  if (!svg) return;
  const clone = svg.cloneNode(true) as SVGSVGElement; clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  downloadText(name, new XMLSerializer().serializeToString(clone), "image/svg+xml;charset=utf-8");
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

function AlignmentCanvas({ sequences, alphabet, rowHeight, cellWidth, selected, onSelect, scrollRef, peerRef }: {
  sequences: string[]; alphabet: "nt" | "aa"; rowHeight: number; cellWidth: number; selected: number; onSelect(index: number): void;
  scrollRef: React.RefObject<HTMLDivElement | null>; peerRef: React.RefObject<HTMLDivElement | null>;
}) {
  const canvas = useRef<HTMLCanvasElement>(null), size = useViewportSize(scrollRef), columns = sequences[0]?.length ?? 0;
  useEffect(() => {
    const container = scrollRef.current, target = canvas.current; if (!container || !target) return;
    const draw = () => {
      const context = canvasResolution(target, size.width, size.height), left = container.scrollLeft, top = container.scrollTop;
      const firstRow = Math.max(0, Math.floor(top / rowHeight)), lastRow = Math.min(sequences.length, Math.ceil((top + size.height) / rowHeight) + 1);
      const firstColumn = Math.max(0, Math.floor(left / cellWidth)), lastColumn = Math.min(columns, Math.ceil((left + size.width) / cellWidth) + 1);
      context.clearRect(0, 0, size.width, size.height); context.textAlign = "center"; context.textBaseline = "middle";
      context.font = `700 ${Math.min(11, rowHeight - 4, cellWidth - 2)}px ui-monospace,monospace`;
      for (let row = firstRow; row < lastRow; row += 1) for (let column = firstColumn; column < lastColumn; column += 1) {
        const base = sequences[row]?.[column] ?? "-", x = column * cellWidth - left, y = row * rowHeight - top;
        context.fillStyle = alphabet === "nt" ? NT_COLORS[base.toUpperCase()] ?? "#c7cbc7" : aaColor(base);
        context.fillRect(x, y, cellWidth + 0.4, rowHeight + 0.4);
        if (cellWidth >= 9) { context.fillStyle = "#132321"; context.fillText(base, x + cellWidth / 2, y + rowHeight / 2); }
      }
      if (selected >= firstRow && selected < lastRow) {
        context.strokeStyle = "#132321"; context.lineWidth = 2; context.strokeRect(0, selected * rowHeight - top, size.width, rowHeight);
      }
    };
    draw();
    const sync = () => {
      if (peerRef.current && Math.abs(peerRef.current.scrollTop - container.scrollTop) > 1) peerRef.current.scrollTop = container.scrollTop;
      draw();
    };
    container.addEventListener("scroll", sync, { passive: true }); return () => container.removeEventListener("scroll", sync);
  }, [sequences, alphabet, rowHeight, cellWidth, selected, columns, size, scrollRef, peerRef]);
  return <div ref={scrollRef} className="alignment-viewport" onPointerDown={(event) => {
    const index = Math.floor((scrollRef.current!.scrollTop + event.nativeEvent.offsetY) / rowHeight);
    if (index >= 0 && index < sequences.length) onSelect(index);
  }}><div style={{ width: Math.max(size.width, columns * cellWidth), height: Math.max(size.height, sequences.length * rowHeight), position: "relative" }}><canvas ref={canvas} /></div></div>;
}

function niceScale(maximum: number): number {
  if (!(maximum > 0)) return 0;
  const target = maximum / 4, power = 10 ** Math.floor(Math.log10(target)), scaled = target / power;
  return (scaled >= 5 ? 5 : scaled >= 2 ? 2 : 1) * power;
}

function TreeSvg({ root, width, rowHeight, selected, onSelect, scrollRef, peerRef, svgRef }: {
  root: TreeNode; width: number; rowHeight: number; selected: number; onSelect(index: number): void;
  scrollRef: React.RefObject<HTMLDivElement | null>; peerRef: React.RefObject<HTMLDivElement | null>;
  svgRef: React.RefObject<SVGSVGElement | null>;
}) {
  const labelWidth = Math.min(210, Math.max(110, width * 0.34)), branchWidth = Math.max(80, width - labelWidth);
  const layout = useMemo(() => layoutTree(root, branchWidth, rowHeight, "phylogram", 18, rowHeight / 2), [root, branchWidth, rowHeight]);
  const leafNodes = useMemo(() => layout.nodes.filter((node) => !node.children.length).sort((left, right) => left.y - right.y), [layout]);
  const scale = niceScale(layout.maximumDistance), horizontalPadding = Math.min(24, branchWidth * 0.1);
  const scalePixels = scale && layout.maximumDistance ? scale / layout.maximumDistance * Math.max(0, branchWidth - 2 * horizontalPadding) : 0;
  useEffect(() => {
    const container = scrollRef.current; if (!container) return;
    const sync = () => { if (peerRef.current && Math.abs(peerRef.current.scrollTop - container.scrollTop) > 1) peerRef.current.scrollTop = container.scrollTop; };
    container.addEventListener("scroll", sync, { passive: true }); return () => container.removeEventListener("scroll", sync);
  }, [scrollRef, peerRef]);
  const height = Math.max(layout.height, leafNodes.length * rowHeight);
  return <div ref={scrollRef} className="tree-viewport" onPointerDown={(event) => {
    const index = Math.floor((scrollRef.current!.scrollTop + event.nativeEvent.offsetY) / rowHeight);
    if (index >= 0 && index < leafNodes.length) onSelect(index);
  }}><svg ref={svgRef} width={width} height={height} viewBox={`0 0 ${width} ${height}`} role="img" aria-label="FastTree phylogram coordinated with aligned sequences">
    <rect width={width} height={height} fill="#fbfaf5" />
    {selected >= 0 && <rect x="0" y={selected * rowHeight} width={width} height={rowHeight} fill="rgba(201,239,115,.32)" />}
    {scalePixels >= 8 && <g><line x1={horizontalPadding} x2={horizontalPadding + scalePixels} y1="9" y2="9" stroke="#273d37" strokeWidth="1.2" /><line x1={horizontalPadding} x2={horizontalPadding} y1="5" y2="13" stroke="#273d37" /><line x1={horizontalPadding + scalePixels} x2={horizontalPadding + scalePixels} y1="5" y2="13" stroke="#273d37" /><text x={horizontalPadding + scalePixels / 2} y="22" textAnchor="middle" fontFamily="Inter,Arial,sans-serif" fontSize="8" fill="#4e5d58">{scale.toPrecision(2)} substitutions/site</text></g>}
    <line x1={branchWidth} x2={branchWidth} y1="0" y2={height} stroke="#d2d8d3" />
    {layout.edges.map((edge, index) => <path key={`edge-${index}`} d={`M${edge.parent.x},${edge.parent.y} V${edge.child.y} H${edge.child.x}`} stroke="#3f5650" strokeWidth="0.95" strokeLinecap="square" strokeLinejoin="miter" fill="none" />)}
    {layout.nodes.map((node, index) => {
      const leaf = !node.children.length;
      return <g key={`node-${index}`}><circle cx={node.x} cy={node.y} r={leaf ? 2.6 : 1.6} fill={leaf ? "#08796f" : "#fbfaf5"} stroke="#31534b" strokeWidth=".75"><title>{node.name || "internal node"}</title></circle>{leaf && <><line x1={node.x + 4} x2={branchWidth + 4} y1={node.y} y2={node.y} stroke="#b3bdb8" strokeWidth=".65" strokeDasharray="2 3" /><text x={branchWidth + 8} y={node.y + 3.3} fontFamily="ui-monospace,monospace" fontSize="9" fill="#263630">{shortName(node.name)}<title>{node.name}</title></text></>}</g>;
    })}
  </svg></div>;
}

export function AlignmentTreeViewer({ fasta, newick, alphabet = "nt", name = "webporpid-tree" }: { fasta: string; newick?: string; alphabet?: "nt" | "aa"; name?: string }) {
  const records = useMemo(() => parseFasta(fasta), [fasta]);
  const [cellWidth, setCellWidth] = useState(11), [rowHeight, setRowHeight] = useState(20), [treeWidth, setTreeWidth] = useState(560);
  const [selected, setSelected] = useState(0), [layoutMode, setLayoutMode] = useState<TreeLayoutMode>("phylogram"), [ladderization, setLadderization] = useState<"none" | LadderizeDirection>("none");
  const treeScroll = useRef<HTMLDivElement>(null), alignmentScroll = useRef<HTMLDivElement>(null), svgRef = useRef<SVGSVGElement>(null);
  const parsed = useMemo(() => {
    if (!newick) return { tree: undefined, error: "" };
    try { return { tree: parseNewick(newick), error: "" }; }
    catch (cause) { return { tree: undefined, error: cause instanceof Error ? cause.message : String(cause) }; }
  }, [newick]);
  const orderedTree = useMemo(() => !parsed.tree || ladderization === "none" ? parsed.tree : ladderizeTree(parsed.tree, ladderization), [parsed.tree, ladderization]);
  const displayTree = useMemo(() => {
    if (!orderedTree || layoutMode === "phylogram") return orderedTree;
    const unit = (node: TreeNode): TreeNode => ({ ...node, length: node.length, children: node.children.map((child) => ({ ...unit(child), length: 1 })) });
    return unit(orderedTree);
  }, [orderedTree, layoutMode]);
  const orderedRecords = useMemo(() => {
    if (!displayTree) return records;
    const byName = new Map(records.map((record) => [safeTreeName(record.name), record]));
    const found = leaves(displayTree).map((leaf) => byName.get(leaf)).filter((record): record is (typeof records)[number] => Boolean(record));
    return found.length === records.length ? found : records;
  }, [displayTree, records]);
  useEffect(() => { setSelected(0); }, [fasta, newick, ladderization]);
  const sequences = orderedRecords.map((record) => record.sequence);
  return <section className="alignment-tree-viewer">
    {parsed.error && <div className="viewer-status">The stored tree could not be rendered ({parsed.error}); the alignment remains available.</div>}
    <header><div><span className="section-kicker">Swig phylogram + virtual alignment</span><h3>{alphabet === "nt" ? "Nucleotide" : "Amino-acid"} variants</h3></div>
      <div className="viewer-controls">
        {displayTree && <><label>Tree width <input type="range" min="260" max="1100" step="10" value={treeWidth} onChange={(event) => setTreeWidth(Number(event.target.value))} /></label><label>Tree <select value={layoutMode} onChange={(event) => setLayoutMode(event.target.value as TreeLayoutMode)}><option value="phylogram">Branch lengths</option><option value="cladogram">Topology</option></select></label><label>Tip order <select value={ladderization} onChange={(event) => setLadderization(event.target.value as typeof ladderization)}><option value="none">Newick order</option><option value="large-first">Large clades first</option><option value="small-first">Small clades first</option></select></label></>}
        <label>Columns <input type="range" min="3" max="18" value={cellWidth} onChange={(event) => setCellWidth(Number(event.target.value))} /></label><label>Rows <input type="range" min="14" max="34" value={rowHeight} onChange={(event) => setRowHeight(Number(event.target.value))} /></label>
        {displayTree && <><button type="button" onClick={() => saveSvg(svgRef.current, `${name}.svg`)}>Tree SVG</button><button type="button" onClick={() => downloadText(`${name}.nwk`, `${serializeNewick(displayTree)};\n`, "text/plain;charset=utf-8")}>Newick</button></>}
      </div></header>
    <div className={`viewer-grid ${displayTree ? "with-tree" : "alignment-only"}`} style={displayTree ? { gridTemplateColumns: `${treeWidth}px minmax(0,1fr)` } : undefined}>
      {displayTree && <TreeSvg root={displayTree} width={treeWidth} rowHeight={rowHeight} selected={selected} onSelect={setSelected} scrollRef={treeScroll} peerRef={alignmentScroll} svgRef={svgRef} />}
      <AlignmentCanvas sequences={sequences} alphabet={alphabet} rowHeight={rowHeight} cellWidth={cellWidth} selected={selected} onSelect={setSelected} scrollRef={alignmentScroll} peerRef={treeScroll} />
    </div>
    <footer><strong>{orderedRecords[selected]?.name ?? "No sequence"}</strong><span>{orderedRecords.length.toLocaleString()} rows × {(sequences[0]?.length ?? 0).toLocaleString()} columns · alignment cells are canvas-virtualized</span></footer>
  </section>;
}
