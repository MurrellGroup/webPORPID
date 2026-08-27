import { useEffect, useMemo, useRef, useState } from "react";
import { parseFasta } from "../config";
import { layoutTree, parseNewick, type LayoutNode } from "../tree";

const NT_COLORS: Record<string, string> = { A: "#78c679", C: "#6baed6", G: "#f2c14e", T: "#ef7668", U: "#ef7668", N: "#c7cbc7", "-": "#ece9df" };
const AA_GROUPS: Array<[string, string]> = [["AILMFWV", "#79b8a4"], ["KRH", "#6baed6"], ["DE", "#ef7668"], ["STNQ", "#c9df75"], ["CGP", "#e5b85c"], ["Y", "#a98dcc"], ["X*-", "#c7cbc7"]];
const aaColor = (base: string) => AA_GROUPS.find(([letters]) => letters.includes(base.toUpperCase()))?.[1] ?? "#ddd8cc";

function useCanvasSize(ref: React.RefObject<HTMLDivElement | null>) {
  const [size, setSize] = useState({ width: 1, height: 1 });
  useEffect(() => {
    const node = ref.current; if (!node) return;
    const update = () => setSize({ width: Math.max(1, node.clientWidth), height: Math.max(1, node.clientHeight) }); update();
    const observer = new ResizeObserver(update); observer.observe(node); return () => observer.disconnect();
  }, [ref]); return size;
}

function canvasResolution(canvas: HTMLCanvasElement, width: number, height: number) {
  const ratio = window.devicePixelRatio || 1; canvas.width = Math.floor(width * ratio); canvas.height = Math.floor(height * ratio);
  canvas.style.width = `${width}px`; canvas.style.height = `${height}px`; const context = canvas.getContext("2d")!; context.setTransform(ratio, 0, 0, ratio, 0, 0); return context;
}

function TreeCanvas({ newick, names, rowHeight, selected, onSelect, scrollRef, peerRef }: {
  newick: string; names: string[]; rowHeight: number; selected: number; onSelect(index: number): void;
  scrollRef: React.RefObject<HTMLDivElement | null>; peerRef: React.RefObject<HTMLDivElement | null>;
}) {
  const canvas = useRef<HTMLCanvasElement>(null), size = useCanvasSize(scrollRef);
  const safeNames = useMemo(() => names.map((name) => name.replace(/[^A-Za-z0-9_.|*+\-]/g, "_")), [names]);
  const rows = useMemo(() => new Map(safeNames.map((name, index) => [name, index])), [safeNames]);
  const tree = useMemo(() => layoutTree(parseNewick(newick), rows, size.width, rowHeight), [newick, rows, size.width, rowHeight]);
  useEffect(() => {
    const container = scrollRef.current, target = canvas.current; if (!container || !target) return;
    const draw = () => {
      const context = canvasResolution(target, size.width, size.height), top = container.scrollTop, bottom = top + size.height;
      context.clearRect(0, 0, size.width, size.height); context.translate(0, -top); context.strokeStyle = "#71817b"; context.lineWidth = 1;
      const visit = (node: LayoutNode) => {
        if (node.children.length) {
          const ys = node.children.map((child) => child.y); context.beginPath(); context.moveTo(node.x, Math.min(...ys)); context.lineTo(node.x, Math.max(...ys)); context.stroke();
        }
        for (const child of node.children) {
          if (child.y >= top - rowHeight && child.y <= bottom + rowHeight) { context.beginPath(); context.moveTo(node.x, child.y); context.lineTo(child.x, child.y); context.stroke(); }
          visit(child);
        }
        if (!node.children.length && node.y >= top - rowHeight && node.y <= bottom + rowHeight) {
          const index = rows.get(node.name.replace(/[^A-Za-z0-9_.|*+\-]/g, "_")) ?? -1;
          if (index === selected) { context.fillStyle = "rgba(201,239,115,.36)"; context.fillRect(0, node.y - rowHeight / 2, size.width, rowHeight); }
          context.fillStyle = "#132321"; context.font = "11px ui-monospace, monospace"; context.textBaseline = "middle";
          context.fillText(names[index] ?? node.name, Math.min(node.x + 5, size.width - 112), node.y, 108);
        }
      }; visit(tree);
    };
    draw(); const sync = () => { if (peerRef.current && Math.abs(peerRef.current.scrollTop - container.scrollTop) > 1) peerRef.current.scrollTop = container.scrollTop; draw(); };
    container.addEventListener("scroll", sync, { passive: true }); return () => container.removeEventListener("scroll", sync);
  }, [tree, rows, names, rowHeight, selected, size, scrollRef, peerRef]);
  return <div ref={scrollRef} className="tree-viewport" onPointerDown={(event) => {
    const index = Math.floor((scrollRef.current!.scrollTop + event.nativeEvent.offsetY) / rowHeight); if (index >= 0 && index < names.length) onSelect(index);
  }}><div style={{ height: names.length * rowHeight, minHeight: size.height, position: "relative" }}><canvas ref={canvas} /></div></div>;
}

function AlignmentCanvas({ sequences, alphabet, rowHeight, cellWidth, selected, onSelect, scrollRef, peerRef }: {
  sequences: string[]; alphabet: "nt" | "aa"; rowHeight: number; cellWidth: number; selected: number; onSelect(index: number): void;
  scrollRef: React.RefObject<HTMLDivElement | null>; peerRef: React.RefObject<HTMLDivElement | null>;
}) {
  const canvas = useRef<HTMLCanvasElement>(null), size = useCanvasSize(scrollRef), columns = sequences[0]?.length ?? 0;
  useEffect(() => {
    const container = scrollRef.current, target = canvas.current; if (!container || !target) return;
    const draw = () => {
      const context = canvasResolution(target, size.width, size.height), left = container.scrollLeft, top = container.scrollTop;
      const firstRow = Math.max(0, Math.floor(top / rowHeight)), lastRow = Math.min(sequences.length, Math.ceil((top + size.height) / rowHeight));
      const firstColumn = Math.max(0, Math.floor(left / cellWidth)), lastColumn = Math.min(columns, Math.ceil((left + size.width) / cellWidth));
      context.clearRect(0, 0, size.width, size.height); context.textAlign = "center"; context.textBaseline = "middle"; context.font = `${Math.min(12, rowHeight - 3)}px ui-monospace, monospace`;
      for (let row = firstRow; row < lastRow; row++) for (let column = firstColumn; column < lastColumn; column++) {
        const base = sequences[row][column] ?? "-", x = column * cellWidth - left, y = row * rowHeight - top;
        context.fillStyle = alphabet === "nt" ? NT_COLORS[base.toUpperCase()] ?? "#c7cbc7" : aaColor(base); context.fillRect(x, y, cellWidth + 0.4, rowHeight + 0.4);
        if (cellWidth >= 9) { context.fillStyle = "#132321"; context.fillText(base, x + cellWidth / 2, y + rowHeight / 2); }
      }
      if (selected >= firstRow && selected < lastRow) { context.strokeStyle = "#132321"; context.lineWidth = 2; context.strokeRect(0, selected * rowHeight - top, size.width, rowHeight); }
    };
    draw(); const sync = () => { if (peerRef.current && Math.abs(peerRef.current.scrollTop - container.scrollTop) > 1) peerRef.current.scrollTop = container.scrollTop; draw(); };
    container.addEventListener("scroll", sync, { passive: true }); return () => container.removeEventListener("scroll", sync);
  }, [sequences, alphabet, rowHeight, cellWidth, selected, columns, size, scrollRef, peerRef]);
  return <div ref={scrollRef} className="alignment-viewport" onPointerDown={(event) => {
    const index = Math.floor((scrollRef.current!.scrollTop + event.nativeEvent.offsetY) / rowHeight); if (index >= 0 && index < sequences.length) onSelect(index);
  }}><div style={{ width: Math.max(size.width, columns * cellWidth), height: Math.max(size.height, sequences.length * rowHeight), position: "relative" }}><canvas ref={canvas} /></div></div>;
}

export function AlignmentTreeViewer({ fasta, newick, alphabet = "nt" }: { fasta: string; newick?: string; alphabet?: "nt" | "aa" }) {
  const rows = useMemo(() => parseFasta(fasta), [fasta]), [cellWidth, setCellWidth] = useState(11), [rowHeight, setRowHeight] = useState(20), [selected, setSelected] = useState(0);
  const treeScroll = useRef<HTMLDivElement>(null), alignmentScroll = useRef<HTMLDivElement>(null), sequences = rows.map((row) => row.sequence);
  return <section className="alignment-tree-viewer">
    <header><div><span className="section-kicker">Linked phylogeny + alignment</span><h3>{alphabet === "nt" ? "Nucleotide" : "Amino-acid"} variants</h3></div>
      <div className="viewer-controls"><label>Columns <input type="range" min="3" max="18" value={cellWidth} onChange={(event) => setCellWidth(Number(event.target.value))} /></label><label>Rows <input type="range" min="14" max="30" value={rowHeight} onChange={(event) => setRowHeight(Number(event.target.value))} /></label></div></header>
    <div className={`viewer-grid ${newick ? "with-tree" : "alignment-only"}`}>
      {newick && <TreeCanvas newick={newick} names={rows.map((row) => row.name)} rowHeight={rowHeight} selected={selected} onSelect={setSelected} scrollRef={treeScroll} peerRef={alignmentScroll} />}
      <AlignmentCanvas sequences={sequences} alphabet={alphabet} rowHeight={rowHeight} cellWidth={cellWidth} selected={selected} onSelect={setSelected} scrollRef={alignmentScroll} peerRef={treeScroll} />
    </div>
    <footer><strong>{rows[selected]?.name ?? "No sequence"}</strong><span>{rows.length.toLocaleString()} rows × {(rows[0]?.sequence.length ?? 0).toLocaleString()} columns · canvas-virtualized</span></footer>
  </section>;
}
