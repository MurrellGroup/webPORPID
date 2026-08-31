import { inputFilterStats, parameterSettings, sampleOverviewStats } from "./report-stats.ts";
import { OPTIONAL_STAGE_ORDER, stageCompleted, stageStatus } from "./optional-stages.ts";
import type { ResultBundle } from "./types.ts";

const quote = (value: unknown) => {
  const source = value == null ? "" : String(value);
  return /[",\r\n]/.test(source) ? `"${source.replaceAll('"', '""')}"` : source;
};
const csv = (headers: string[], rows: unknown[][]) => [headers, ...rows].map((row) => row.map(quote).join(",")).join("\n") + "\n";

export interface OverviewExport { name: string; text: string }

function thresholdDecisionRows(bundle: ResultBundle): unknown[][] {
  const rows: unknown[][] = [];
  for (const record of bundle.thresholdSelections ?? []) {
    for (const [parameter, value] of Object.entries(record.parameters))
      rows.push([record.id, record.phase, record.acceptedUtc, "global", "", parameter, value, ""]);
    for (const sample of record.samples) {
      const parameters = record.phase === "umi"
        ? [["familySizeOverride", sample.familySizeOverride]] as const
        : [["artefactFractionOverride", sample.artefactFractionOverride], ["outlierQuantileOverride", sample.outlierQuantileOverride],
          ["agreementOverride", sample.agreementOverride]] as const;
      for (const [parameter, value] of parameters)
        rows.push([record.id, record.phase, record.acceptedUtc, "sample", sample.sample, parameter, value ?? "global", ""]);
    }
    for (const change of record.changes)
      rows.push([record.id, record.phase, record.acceptedUtc, "audit", "", "", "", change]);
  }
  return rows;
}

/** Machine-readable counterparts of every run-wide overview section. */
export function buildOverviewExports(bundle: ResultBundle): OverviewExport[] {
  const overview = sampleOverviewStats(bundle), input = inputFilterStats(bundle), parameters = parameterSettings(bundle);
  const contaminationComplete = stageCompleted(bundle, "contamination"), postprocessingComplete = stageCompleted(bundle, "postprocessing"), collapseComplete = stageCompleted(bundle, "collapse");
  const postprocessingFields = new Set(["retainedFamilies", "functionalEvaluatedFamilies", "functionalPassedFamilies", "artefactReadPercent",
    "agreementReadPercent", "panelReadPercent", "functionalReadPercent"]);
  const summaryKeys = ["sample", "donorId", "demultiplexedReads", "selectedReads", "downsampledReads", "downsampledPercent", "observedFamilies",
    "consensusFamilies", "retainedFamilies", "functionalConfigured", "functionalEvaluatedFamilies", "functionalPassedFamilies",
    "collapsedHaplotypes", "bpbReadPercent", "umiLengthReadPercent", "familySizeReadPercent", "ldaReadPercent",
    "heteroduplexReadPercent", "artefactReadPercent", "agreementReadPercent", "contaminationReadPercent", "panelReadPercent",
    "functionalReadPercent"] as const;
  return [
    { name: "README.txt", text: "Cross-sample tables from the webPORPID overview. Percentage columns are percentages (0–100), not fractions. Optional-stage status distinguishes zero results from work that was not computed.\n" },
    { name: "sample-summary.csv", text: csv([...summaryKeys], overview.map((row) => summaryKeys.map((key) =>
      key === "contaminationReadPercent" && !contaminationComplete ? undefined
        : key === "collapsedHaplotypes" && !collapseComplete ? undefined
          : postprocessingFields.has(key) && !postprocessingComplete ? undefined : row[key]))) },
    { name: "input-filtering.csv", text: csv(["key", "label", "count", "percent", "note"], input.map((row) => [row.key, row.label, row.count, row.percent, row.note])) },
    { name: "parameters.csv", text: csv(["scope", "sample", "parameter", "value"], parameters.map((row) => [row.scope, row.sample, row.parameter, row.value])) },
    { name: "optional-stage-status.csv", text: csv(["stage", "state", "detail", "updated_utc"], OPTIONAL_STAGE_ORDER.map((stage) => {
      const status = stageStatus(bundle, stage); return [stage, status.state, status.detail, status.updatedUtc];
    })) },
    { name: "timings.csv", text: csv(["stage", "seconds", "work_items"], (bundle.timings ?? []).map((row) => [row.stage, row.seconds, row.workItems])) },
    { name: "interactive-threshold-decisions.csv", text: csv(
      ["checkpoint_id", "phase", "accepted_utc", "scope", "sample", "parameter", "value", "change"], thresholdDecisionRows(bundle)) },
    { name: "input-file-mappings.csv", text: csv(["slot", "role", "expected_name", "uploaded_name", "uploaded_size"],
      (bundle.inputMappings ?? []).map((row) => [row.slot, row.role, row.expectedName, row.uploadedName, row.uploadedSize])) },
    { name: "provenance.csv", text: csv(["field", "value"], Object.entries(bundle.provenance).map(([key, value]) => [key, value])) },
  ];
}
