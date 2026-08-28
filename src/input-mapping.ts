import type { InputFileMapping, PipelineConfig } from "./types";

export interface ReferenceSlot {
  id: string;
  expectedName: string;
  role: InputFileMapping["role"];
  label: string;
  samples: string[];
}

export const normalizeInputPath = (value: string) => value.replaceAll("\\", "/").replace(/^\.\//, "");
export const inputBasename = (value: string) => normalizeInputPath(value).split("/").at(-1)!.toLowerCase();

export function referenceSlots(config: PipelineConfig): ReferenceSlot[] {
  const slots = new Map<string, ReferenceSlot>();
  const add = (expectedName: string, role: ReferenceSlot["role"], label: string, sample?: string) => {
    const id = normalizeInputPath(expectedName);
    const existing = slots.get(id);
    if (existing) {
      if (sample && !existing.samples.includes(sample)) existing.samples.push(sample);
      return;
    }
    slots.set(id, { id, expectedName, role, label, samples: sample ? [sample] : [] });
  };
  for (const sample of config.samples) {
    add(sample.panel, "panel", "Panel", sample.name);
    if (sample.functionalReference) add(sample.functionalReference, "functional-reference", "Functional reference", sample.name);
  }
  if (config.parameters.contaminationFilter) add(config.contaminationPanel, "contamination-panel", "Contamination panel");
  return [...slots.values()];
}

export function nameMatchingSlot(fileName: string, slots: readonly ReferenceSlot[], excludedId?: string): ReferenceSlot | undefined {
  const base = inputBasename(fileName);
  const matches = slots.filter((slot) => slot.id !== excludedId && inputBasename(slot.expectedName) === base);
  return matches.length === 1 ? matches[0] : undefined;
}

export function referenceFileMap(assignments: Readonly<Record<string, File>>): Map<string, () => Promise<string>> {
  return new Map(Object.entries(assignments).map(([slot, file]) => [slot, () => file.text()]));
}

export function referenceMappingRecords(slots: readonly ReferenceSlot[], assignments: Readonly<Record<string, File>>): InputFileMapping[] {
  return slots.flatMap((slot) => {
    const file = assignments[slot.id];
    return file ? [{ slot: slot.id, role: slot.role, expectedName: slot.expectedName, uploadedName: file.name, uploadedSize: file.size }] : [];
  });
}
