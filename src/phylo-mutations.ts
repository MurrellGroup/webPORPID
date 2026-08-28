import { translateAlignedNucleotides } from "./alignment-utils.ts";
import type { TreeNode } from "./tree.ts";

export interface BranchMutation {
  column: number;
  from: string;
  to: string;
  childClade: string;
}

export interface ParsimonyMap {
  score: number;
  sequencesByClade: Map<string, string>;
  mutationsByClade: Map<string, BranchMutation[]>;
}

const STATES = ["A", "C", "G", "T", "-"] as const;
const MASKS: Record<string, number> = {
  A: 1, C: 2, G: 4, T: 8, U: 8, R: 5, Y: 10, S: 6, W: 9, K: 12, M: 3,
  B: 14, D: 13, H: 11, V: 7, N: 15, X: 15, "?": 15, "-": 16, ".": 16,
};
const INF = 1_000_000;

function leafNames(node: TreeNode, target: string[] = []): string[] {
  if (!node.children.length) target.push(node.name);
  else node.children.forEach((child) => leafNames(child, target));
  return target;
}

/** Stable node identity independent of child order. */
export function cladeSignature(node: TreeNode): string {
  return leafNames(node).sort().join("\u0000");
}

function mask(value: string): number { return MASKS[value.toUpperCase()] ?? 15; }
function choose(costs: readonly number[], preferred = -1): number {
  const best = Math.min(...costs);
  if (preferred >= 0 && costs[preferred] === best) return preferred;
  return costs.findIndex((value) => value === best);
}

/** Equal-cost Sankoff reconstruction with gaps represented as an explicit state. */
export function mapParsimonyMutations(root: TreeNode, alignedByName: ReadonlyMap<string, string>): ParsimonyMap {
  if (!alignedByName.size) throw new Error("Mutation mapping requires an alignment.");
  const widths = new Set([...alignedByName.values()].map((sequence) => sequence.length));
  if (widths.size !== 1) throw new Error("Mutation mapping requires equal-length sequences.");
  const columns = widths.values().next().value as number;
  const nodes: TreeNode[] = [];
  const collect = (node: TreeNode) => { nodes.push(node); node.children.forEach(collect); };
  collect(root);
  const sequenceBuffers = new Map(nodes.map((node) => [node, [] as string[]]));
  const mutationBuffers = new Map<string, BranchMutation[]>();
  let score = 0;

  for (let column = 0; column < columns; column += 1) {
    const costs = new Map<TreeNode, number[]>();
    const postorder = (node: TreeNode): number[] => {
      if (!node.children.length) {
        const observed = mask(alignedByName.get(node.name)?.[column] ?? "N");
        const result = STATES.map((_, state) => observed & (1 << state) ? 0 : INF);
        costs.set(node, result); return result;
      }
      const childCosts = node.children.map(postorder);
      const result = STATES.map((_, parentState) => childCosts.reduce((sum, child) => {
        let best = INF;
        for (let childState = 0; childState < STATES.length; childState += 1)
          best = Math.min(best, child[childState] + Number(childState !== parentState));
        return sum + best;
      }, 0));
      costs.set(node, result); return result;
    };
    const rootCosts = postorder(root), rootState = choose(rootCosts); score += rootCosts[rootState];
    const preorder = (node: TreeNode, state: number) => {
      sequenceBuffers.get(node)!.push(STATES[state]);
      for (const child of node.children) {
        const childCosts = costs.get(child)!;
        const transitioned = childCosts.map((value, childState) => value + Number(childState !== state));
        const childState = choose(transitioned, state);
        if (childState !== state) {
          const signature = cladeSignature(child), mutations = mutationBuffers.get(signature) ?? [];
          mutations.push({ column, from: STATES[state], to: STATES[childState], childClade: signature });
          mutationBuffers.set(signature, mutations);
        }
        preorder(child, childState);
      }
    };
    preorder(root, rootState);
  }
  const sequencesByClade = new Map<string, string>();
  for (const [node, sequence] of sequenceBuffers) sequencesByClade.set(cladeSignature(node), sequence.join(""));
  return { score, sequencesByClade, mutationsByClade: mutationBuffers };
}

export function aminoAcidBranchMutations(parent: string, child: string, childClade: string, frameOffset = 0): BranchMutation[] {
  const left = translateAlignedNucleotides(parent, frameOffset), right = translateAlignedNucleotides(child, frameOffset), result: BranchMutation[] = [];
  for (let column = 0; column < Math.min(left.length, right.length); column += 1) {
    const from = left[column], to = right[column];
    if (from !== to && from !== "X" && to !== "X") result.push({ column, from, to, childClade });
  }
  return result;
}
