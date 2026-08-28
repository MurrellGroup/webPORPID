export interface TreeNode {
  name: string;
  length: number;
  children: TreeNode[];
}

export interface TreeLayoutNode extends TreeNode {
  x: number;
  y: number;
}

export interface TreeLayout {
  nodes: TreeLayoutNode[];
  edges: Array<{ parent: TreeLayoutNode; child: TreeLayoutNode }>;
  width: number;
  height: number;
  leaves: number;
  mode: TreeLayoutMode;
  maximumDistance: number;
}

export type TreeLayoutMode = "phylogram" | "cladogram";
export type LadderizeDirection = "large-first" | "small-first";

export const FASTTREE_DOUBLE_MINIMUM_BRANCH = 5e-9;
export const FASTTREE_AMBIGUOUS_BRANCH_THRESHOLD = 1e-8;

export interface CollapsedTree {
  root: TreeNode;
  collapsedEdges: number;
  threshold: number;
}

export function extractNewick(output: string): string {
  const text = String(output).replace(/\x1b\[[0-9;]*m/g, "");
  let best = "";
  for (let start = 0; start < text.length; start += 1) {
    if (text[start] !== "(") continue;
    let depth = 0;
    let quote = "";
    let commentDepth = 0;
    for (let index = start; index < text.length; index += 1) {
      const value = text[index];
      if (quote) {
        if (value === quote && text[index - 1] !== "\\") quote = "";
        continue;
      }
      if (value === "'" || value === '"') {
        quote = value;
        continue;
      }
      if (value === "[") {
        commentDepth += 1;
        continue;
      }
      if (value === "]" && commentDepth) {
        commentDepth -= 1;
        continue;
      }
      if (commentDepth) continue;
      if (value === "(") depth += 1;
      else if (value === ")") {
        depth -= 1;
        if (depth !== 0) continue;
        let end = index + 1;
        while (end < text.length && text[end] !== ";" && text[end] !== "\n" && text[end] !== "\r") end += 1;
        if (text[end] !== ";") break;
        end += 1;
        const candidate = text.slice(start, end).trim();
        if ((candidate.includes(",") || candidate.includes(":")) && candidate.length > best.length) best = candidate;
        break;
      }
      if (depth < 0) break;
    }
  }
  if (!best) throw new Error("FastTree did not return a complete Newick tree.");
  return `${best.replace(/;+$/, "")};`;
}

export function parseNewick(text: string): TreeNode {
  const source = extractNewick(text).replace(/\[[^\]]*\]/g, "");
  let position = 0;
  const whitespace = () => {
    while (/\s/.test(source[position] ?? "")) position += 1;
  };
  const token = () => {
    whitespace();
    if (source[position] === "'" || source[position] === '"') {
      const quote = source[position++];
      let value = "";
      while (position < source.length && source[position] !== quote) value += source[position++];
      if (source[position] === quote) position += 1;
      whitespace();
      return value;
    }
    const start = position;
    while (position < source.length && !"(),:;".includes(source[position]) && !/\s/.test(source[position])) position += 1;
    const value = source.slice(start, position).trim();
    whitespace();
    return value;
  };
  const length = () => {
    whitespace();
    if (source[position] !== ":") return 0;
    position += 1;
    const value = token();
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed : 0;
  };
  const parse = (): TreeNode => {
    whitespace();
    const children: TreeNode[] = [];
    if (source[position] === "(") {
      position += 1;
      while (position < source.length) {
        children.push(parse());
        if (source[position] === ",") {
          position += 1;
          continue;
        }
        if (source[position] === ")") {
          position += 1;
          break;
        }
        throw new Error("Could not parse the complete Newick tree.");
      }
    }
    const name = token();
    const branchLength = length();
    return { name, length: branchLength, children };
  };
  const root = parse();
  whitespace();
  if (source[position] === ";") position += 1;
  whitespace();
  if (position < source.length) throw new Error("Could not parse the complete Newick tree.");
  return root;
}

export function serializeNewick(node: TreeNode): string {
  const formatLength = (value: number) => Number.isFinite(value) && value !== 0 ? String(value) : "0";
  const visit = (current: TreeNode, isRoot: boolean): string => {
    const children = current.children.length ? `(${current.children.map((child) => visit(child, false)).join(",")})` : "";
    const name = current.name.replace(/[\s():;,]/g, "_");
    // A missing Newick length means "unspecified", not reliably zero. Emit an
    // explicit value for every edge, including FastTree's zero-length tips.
    const length = isRoot ? "" : `:${formatLength(current.length)}`;
    return `${children}${name}${length}`;
  };
  return visit(node, true);
}

export function rootOnOutgroup(root: TreeNode, outgroupName: string): TreeNode {
  const adjacency = new Map<TreeNode, Array<{ node: TreeNode; length: number }>>();
  let outgroup: TreeNode | null = null;
  const visit = (node: TreeNode, parent?: TreeNode) => {
    if (!node.children.length && node.name === outgroupName) outgroup = node;
    if (!adjacency.has(node)) adjacency.set(node, []);
    for (const child of node.children) {
      adjacency.get(node)!.push({ node: child, length: child.length });
      adjacency.set(child, [{ node, length: child.length }]);
      visit(child, node);
    }
    void parent;
  };
  visit(root);
  if (!outgroup) throw new Error(`The tree does not contain the outgroup ${outgroupName}.`);
  const neighbor = adjacency.get(outgroup)?.[0];
  if (!neighbor) return root;
  const orient = (node: TreeNode, parent: TreeNode, branchLength: number): TreeNode => ({
    name: node.name,
    length: branchLength,
    children: (adjacency.get(node) ?? []).filter((edge) => edge.node !== parent).map((edge) => orient(edge.node, node, edge.length)),
  });
  // The N-masked germline represents the UCA anchor, not an ordinary sampled
  // outgroup. Place the root exactly at that tip: its root edge is zero and
  // the ingroup receives the complete original connecting-edge length.
  return {
    name: "",
    length: 0,
    children: [orient(outgroup, neighbor.node, 0), orient(neighbor.node, outgroup, neighbor.length)],
  };
}

/** Reorder children for display without changing topology or branch lengths. */
export function ladderizeTree(root: TreeNode, direction: LadderizeDirection): TreeNode {
  const visit = (node: TreeNode): { node: TreeNode; leaves: number; key: string } => {
    if (!node.children.length) return { node: { ...node, children: [] }, leaves: 1, key: node.name };
    const children = node.children.map(visit).sort((left, right) => {
      const sizeOrder = direction === "large-first" ? right.leaves - left.leaves : left.leaves - right.leaves;
      return sizeOrder || left.key.localeCompare(right.key);
    });
    return {
      node: { ...node, children: children.map((child) => child.node) },
      leaves: children.reduce((sum, child) => sum + child.leaves, 0),
      key: children.map((child) => child.key).sort().join("\u0000"),
    };
  };
  return visit(root).node;
}

/**
 * Return a deterministic child order without changing any splits or lengths.
 * Newick child order is arbitrary, but a stable order prevents equivalent
 * trees from jumping around in the lineage viewer between runtimes.
 */
export function canonicalizeTree(root: TreeNode): TreeNode {
  const visit = (node: TreeNode): { node: TreeNode; key: string } => {
    if (!node.children.length) return { node: { ...node, children: [] }, key: node.name };
    const children = node.children.map(visit).sort((left, right) => left.key.localeCompare(right.key));
    return {
      node: { ...node, children: children.map((child) => child.node) },
      key: children.map((child) => child.key).sort().join("\0"),
    };
  };
  return visit(root).node;
}

/**
 * Collapse internal branches at FastTreeDbl's numerical floor. FastTree can
 * resolve these zero-information edges differently across native and WASM
 * floating-point targets. Leaves are never collapsed, and the untouched raw
 * FastTree Newick remains available separately for audit.
 */
export function collapseShortInternalBranches(
  root: TreeNode,
  threshold = FASTTREE_AMBIGUOUS_BRANCH_THRESHOLD,
): CollapsedTree {
  if (!Number.isFinite(threshold) || threshold < 0) throw new Error("The branch-collapse threshold must be a non-negative number.");
  let collapsedEdges = 0;
  const visit = (node: TreeNode): TreeNode => {
    const children = node.children.flatMap((child) => {
      const next = visit(child);
      if (next.children.length && next.length <= threshold) {
        collapsedEdges += 1;
        return next.children;
      }
      return [next];
    });
    return { ...node, children };
  };
  return { root: canonicalizeTree(visit(root)), collapsedEdges, threshold };
}

export function layoutTree(
  root: TreeNode,
  width = 900,
  rowHeight = 24,
  requestedMode: TreeLayoutMode = "phylogram",
  rightPadding = 220,
  topPadding = 28,
): TreeLayout {
  let hasPositiveLength = false;
  const detectLengths = (node: TreeNode) => {
    for (const child of node.children) {
      if (child.length > 0) hasPositiveLength = true;
      detectLengths(child);
    }
  };
  detectLengths(root);
  const mode: TreeLayoutMode = requestedMode === "phylogram" && hasPositiveLength ? "phylogram" : "cladogram";
  const distances = new Map<TreeNode, number>();
  let maximum = 0;
  const measure = (node: TreeNode, distance: number) => {
    distances.set(node, distance);
    maximum = Math.max(maximum, distance);
    // A zero-length FastTree edge is biologically meaningful. Never inflate it
    // to one unit in a phylogram; use unit depth only in explicit cladogram mode.
    node.children.forEach((child) => measure(child, distance + (mode === "phylogram" ? Math.max(0, child.length) : 1)));
  };
  measure(root, 0);
  const leaves = [...distances.keys()].filter((node) => !node.children.length);
  const leafY = new Map(leaves.map((leaf, index) => [leaf, topPadding + index * rowHeight]));
  const nodes: TreeLayoutNode[] = [];
  const edges: Array<{ parent: TreeLayoutNode; child: TreeLayoutNode }> = [];
  // Padding must contract with the user-controlled width. Otherwise a small
  // width slider value merely hits an invisible 48 px floor instead of truly
  // compressing the tree.
  const leftPadding = Math.min(24, Math.max(0, width * 0.1));
  const effectiveRightPadding = Math.min(Math.max(0, rightPadding), Math.max(0, width * 0.1));
  const drawableWidth = Math.max(0, width - leftPadding - effectiveRightPadding);
  // Branch lengths are normally fractions of one substitution per site. The
  // previous Math.max(maximum, 1) denominator therefore compressed ordinary
  // lineage trees into only a few pixels at the left of the panel. Normalize
  // by the observed maximum itself so the most distant tip uses the full tree
  // viewport, irrespective of the absolute evolutionary scale.
  const distanceScale = maximum > 0 ? maximum : 1;
  const build = (node: TreeNode): TreeLayoutNode => {
    const children = node.children.map(build);
    const y = children.length ? children.reduce((sum, child) => sum + child.y, 0) / children.length : leafY.get(node) ?? 0;
    const layout: TreeLayoutNode = {
      ...node,
      children,
      x: leftPadding + (distances.get(node) ?? 0) / distanceScale * drawableWidth,
      y,
    };
    nodes.push(layout);
    children.forEach((child) => edges.push({ parent: layout, child }));
    return layout;
  };
  build(root);
  return {
    nodes,
    edges,
    width,
    height: Math.max(90, topPadding + Math.max(0, leaves.length - 1) * rowHeight + 32),
    leaves: leaves.length,
    mode,
    maximumDistance: maximum,
  };
}
