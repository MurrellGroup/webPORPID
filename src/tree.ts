export interface TreeNode { name: string; length: number; children: TreeNode[] }

export function parseNewick(source: string): TreeNode {
  source = source.replace(/\[[^\]]*\]/g, "").trim(); let position = 0;
  const whitespace = () => { while (/\s/.test(source[position] ?? "")) position++; };
  const token = () => {
    whitespace(); if (source[position] === "'" || source[position] === '"') {
      const quote = source[position++]; let value = ""; while (position < source.length && source[position] !== quote) value += source[position++];
      if (source[position] === quote) position++; return value;
    }
    const start = position; while (position < source.length && !"(),:;".includes(source[position]) && !/\s/.test(source[position])) position++;
    return source.slice(start, position).trim();
  };
  const branch = () => { whitespace(); if (source[position] !== ":") return 0; position++; const value = Number(token()); return Number.isFinite(value) ? Math.max(0, value) : 0; };
  const parse = (): TreeNode => {
    whitespace(); const children: TreeNode[] = [];
    if (source[position] === "(") {
      position++; while (position < source.length) { children.push(parse()); if (source[position] === ",") { position++; continue; } if (source[position] === ")") { position++; break; } throw new Error("Invalid Newick topology."); }
    }
    return { name: token(), length: branch(), children };
  };
  const root = parse(); whitespace(); if (source[position] === ";") position++; whitespace();
  if (position !== source.length) throw new Error("Trailing content after Newick tree."); return root;
}

export interface LayoutNode extends TreeNode { x: number; y: number; children: LayoutNode[] }
export function layoutTree(root: TreeNode, tipRows: Map<string, number>, width: number, rowHeight: number): LayoutNode {
  let maximum = 0; const distances = new Map<TreeNode, number>();
  const measure = (node: TreeNode, distance: number) => { distances.set(node, distance); maximum = Math.max(maximum, distance); node.children.forEach((child) => measure(child, distance + child.length)); };
  measure(root, 0); const scale = maximum > 0 ? maximum : 1;
  const visit = (node: TreeNode): LayoutNode => {
    const children = node.children.map(visit); const safe = node.name.replace(/[^A-Za-z0-9_.|*+\-]/g, "_");
    const y = children.length ? children.reduce((sum, child) => sum + child.y, 0) / children.length : ((tipRows.get(safe) ?? 0) + 0.5) * rowHeight;
    return { ...node, children, x: 8 + (distances.get(node)! / scale) * (width - 120), y };
  };
  return visit(root);
}
