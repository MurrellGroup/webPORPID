import { useEffect, useMemo, useRef, useState } from "react";
import {
  ALIVIBE_BRIDGE_VERSION, ALIVIBE_SOURCE_REVISION, assertAlivibeInitialLoad, assertAlivibeRoundTripTarget,
  getAlivibeBridge, loadAlivibeNucleotideFasta, readAlivibeNucleotideFasta, type AlivibeEditorWindow, type AlivibeMsaJob,
} from "../alivibe-roundtrip";
import { runAlivibeMsa } from "../alivibe-msa-runtime";
import { alignmentKey, effectiveAlignment, inspectAlignment, summarizeAlignmentChanges, validateCorrectedAlignment, type AlignmentFrameOffset, type AlignmentVariant } from "../alignment-utils";
import { runFastTree } from "../biowasm";
import { parseFasta } from "../config";
import { deduplicateContaminationCalls } from "../contamination";
import { computeThrough, markOptionalStageSkipped, type DeferredAnalysisProgress } from "../deferred-analysis";
import { buildExportArchive } from "../export-archive";
import { OPTIONAL_STAGE_ORDER, stageCompleted, stageStatus } from "../optional-stages";
import { functionalFilterStats, inputFilterStats, parameterSettings, porpidCallStats, postprocFilterStats, sampleOverviewStats, type CountStat, type DualCountStat, type ParameterSettingRow, type SampleOverviewStat } from "../report-stats";
import { exportComponent, type ExportKind, safeDatasetName } from "../result-file";
import { runScalableMsa } from "../scalable-msa";
import type { AlignmentAuditEntry, AlignmentChangeSummary, CollapseGroup, ContaminationCall, OptionalStageName, PostprocRecord, ResultBundle, UmiFamily } from "../types";
import { AgreementPositionPlot, ArtefactDecisionPlot, DinucleotideHeatmaps, MdsApobecPlot, UmiDecisionPlot } from "./charts";
import { AlignmentTreeViewer, type LeafMetadata } from "./alignment-tree-viewer";
import { MethodLink } from "./method-link";

type Tab = "overview" | "families" | "sequences" | "contamination" | "alignment" | "log";
type SortDirection = "asc" | "desc";
interface SortState<K extends string> { key: K; direction: SortDirection }
const PAGE_SIZE = 250;
const CONTAMINATION_TIP_LEGEND = [
  { label: "Contamination-panel reference", color: "#d49a19" },
  { label: "Donor detected contaminant", color: "#c5534f" },
  { label: "Donor non-contaminant", color: "#08796f" },
] as const;

interface AlivibeSession {
  token: string;
  popup: AlivibeEditorWindow;
  sample: string;
  variant: AlignmentVariant;
  baseline: string;
  baselineFingerprint: string;
  frameOffset: AlignmentFrameOffset;
  timer?: number;
}
interface ContaminationPhylogeny {
  sample: string;
  fasta: string;
  newick: string;
  leafMetadata: Record<string, LeafMetadata>;
}
interface DonorPhylogeny {
  donorId: string;
  variant: "collapsed" | "functional";
  fasta: string;
  newick: string;
  leafMetadata: Record<string, LeafMetadata>;
  legend: Array<{ label: string; color: string }>;
}
function downloadData(name: string, data: string | Uint8Array, mime: string) {
  const body: BlobPart = typeof data === "string" ? data : data.slice().buffer as ArrayBuffer;
  const url = URL.createObjectURL(new Blob([body], { type: mime })), anchor = document.createElement("a");
  anchor.href = url; anchor.download = name; anchor.click(); setTimeout(() => URL.revokeObjectURL(url), 0);
}

function Pager({ page, count, onChange }: { page: number; count: number; onChange(page: number): void }) {
  const pages = Math.max(1, Math.ceil(count / PAGE_SIZE)), current = Math.min(page, pages - 1);
  if (count <= PAGE_SIZE) return null;
  return <div className="table-pager"><span>Rows {(current * PAGE_SIZE + 1).toLocaleString()}–{Math.min(count, (current + 1) * PAGE_SIZE).toLocaleString()} of {count.toLocaleString()}</span><div><button type="button" disabled={current === 0} onClick={() => onChange(current - 1)}>Previous</button><strong>{current + 1} / {pages}</strong><button type="button" disabled={current + 1 >= pages} onClick={() => onChange(current + 1)}>Next</button></div></div>;
}

function comparable(value: unknown): string | number {
  if (value == null) return Number.NEGATIVE_INFINITY;
  return typeof value === "number" ? value : typeof value === "boolean" ? Number(value) : String(value).toLowerCase();
}

function sorted<T>(rows: readonly T[], state: SortState<string>, value: (row: T, key: string) => unknown): T[] {
  const direction = state.direction === "asc" ? 1 : -1;
  return [...rows].sort((left, right) => {
    const a = comparable(value(left, state.key)), b = comparable(value(right, state.key));
    return (typeof a === "number" && typeof b === "number" ? a - b : String(a).localeCompare(String(b), undefined, { numeric: true })) * direction;
  });
}

function SortHeader<K extends string>({ label, column, state, onChange }: { label: string; column: K; state: SortState<K>; onChange(state: SortState<K>): void }) {
  const active = state.key === column;
  return <th aria-sort={!active ? "none" : state.direction === "asc" ? "ascending" : "descending"}><button type="button" onClick={() => onChange({ key: column, direction: active && state.direction === "asc" ? "desc" : "asc" })}>{label}<span>{active ? state.direction === "asc" ? "▲" : "▼" : "↕"}</span></button></th>;
}

function Filters({ bundle, sample }: { bundle: ResultBundle; sample: string }) {
  const rows = bundle.records.filter((row) => row.sample === sample), total = rows.length;
  const steps = [
    ["Consensus", total], ["Artefact", rows.filter((row) => row.artefactPass).length],
    ["Agreement", rows.filter((row) => row.artefactPass && row.agreementPass).length],
    ["Contamination", rows.filter((row) => row.artefactPass && row.agreementPass && row.contaminationPass).length],
    ["Panel", rows.filter((row) => row.artefactPass && row.agreementPass && row.contaminationPass && row.panelPass).length],
  ] as const;
  return <div className="filter-flow">{steps.map(([label, count], index) => <div key={label}><span>{label}</span><strong>{count.toLocaleString()}</strong><small>{total ? `${(count / total * 100).toFixed(1)}%` : "—"}</small>{index < steps.length - 1 && <i>→</i>}</div>)}</div>;
}

function ExportMenu({ bundle, sample, allOnly = false }: { bundle: ResultBundle; sample: string; allOnly?: boolean }) {
  const [kind, setKind] = useState<ExportKind>("consensus-fasta"), [exportingAll, setExportingAll] = useState(false);
  const labels: Array<[ExportKind, string]> = [["consensus-fasta", "Consensus FASTA"], ["passed-consensus-fasta", "Passed consensus FASTA"],
    ["rejected-consensus-fasta", "Rejected consensus FASTA"], ["trimmed-nt-fasta", "Collapsed functional nucleotide FASTA"], ["trimmed-aa-fasta", "Collapsed functional amino-acid FASTA"],
    ["family-csv", "UMI family CSV"], ["low-agreement-csv", "Low-agreement CSV"], ["contamination-csv", "Contamination CSV"],
    ["postproc-csv", "Postproc CSV"], ["apobec-csv", "APOBEC CSV"], ["collapse-csv", "Collapse membership CSV"],
    ["nucleotide-alignment", "Collapsed nucleotide alignment"], ["protein-alignment", "Collapsed protein alignment"], ["newick", "Collapsed Newick tree"],
    ["uncollapsed-nucleotide-alignment", "Uncollapsed nucleotide alignment"], ["uncollapsed-protein-alignment", "Uncollapsed protein alignment"], ["uncollapsed-newick", "Uncollapsed Newick tree"],
    ["functional-nucleotide-alignment", "Functional nucleotide alignment"], ["functional-protein-alignment", "Functional protein alignment"],
    ["functional-newick", "Functional Newick tree"], ["log", "Run log"]];
  return <div className="export-menu">{!allOnly && <><select value={kind} onChange={(event) => setKind(event.target.value as ExportKind)}>{labels.map(([value, label]) => <option value={value} key={value}>{label}</option>)}</select>
    <button type="button" onClick={() => { const result = exportComponent(bundle, kind, sample); downloadData(`${safeDatasetName(bundle.config.dataset)}-${safeDatasetName(sample)}.${result.extension}`, result.text, result.mime); }}>Export</button></>}
    <button type="button" disabled={exportingAll} onClick={() => {
      setExportingAll(true);
      window.requestAnimationFrame(() => {
        try { downloadData(`${safeDatasetName(bundle.config.dataset)}-all-outputs.tar.gz`, buildExportArchive(bundle), "application/gzip"); }
        catch (cause) { window.alert(`The export bundle could not be created. ${cause instanceof Error ? cause.message : String(cause)}`); }
        finally { setExportingAll(false); }
      });
    }}>{exportingAll ? "Bundling…" : "Export all (.tar.gz)"}</button></div>;
}

function TimingSummary({ bundle }: { bundle: ResultBundle }) {
  if (!bundle.timings?.length) return null;
  return <article className="timing-card"><header><h3>Run timing</h3><p>Wall time measured around each stored pipeline stage.</p></header><div>{bundle.timings.map((entry) => <span key={entry.stage}><strong>{entry.stage.replaceAll("-", " ")}</strong><em>{entry.seconds < 1 ? `${(entry.seconds * 1000).toFixed(1)} ms` : `${entry.seconds.toFixed(2)} s`}</em></span>)}</div></article>;
}

const pct = (value: number) => `${value.toFixed(1)}%`;

function DualStatsTable({ title, description, rows, familyLabel = "UMI families", readLabel = "Reads / CCS" }: {
  title: string; description: string; rows: DualCountStat[]; familyLabel?: string; readLabel?: string;
}) {
  const [sort, setSort] = useState<SortState<string>>({ key: "label", direction: "asc" });
  const displayed = useMemo(() => sorted(rows, sort, (row, key) => row[key as keyof DualCountStat]), [rows, sort]);
  return <article className="statistics-card"><header><h3>{title}</h3><p>{description}</p></header><div className="table-scroll compact-table"><table><thead><tr><SortHeader label="Call / filter result" column="label" state={sort} onChange={setSort} /><SortHeader label={`${familyLabel} · n`} column="families" state={sort} onChange={setSort} /><SortHeader label={`${familyLabel} · %`} column="familyPercent" state={sort} onChange={setSort} /><SortHeader label={`${readLabel} · n`} column="reads" state={sort} onChange={setSort} /><SortHeader label={`${readLabel} · %`} column="readPercent" state={sort} onChange={setSort} /></tr></thead><tbody>{displayed.map((row) => <tr key={row.key}><td>{row.label}{row.note && <small className="cell-note">{row.note}</small>}</td><td>{row.families.toLocaleString()}</td><td>{pct(row.familyPercent)}</td><td>{row.reads.toLocaleString()}</td><td>{pct(row.readPercent)}</td></tr>)}</tbody></table></div></article>;
}

function CountStatsTable({ rows }: { rows: CountStat[] }) {
  const [sort, setSort] = useState<SortState<string>>({ key: "label", direction: "asc" });
  const displayed = useMemo(() => sorted(rows, sort, (row, key) => row[key as keyof CountStat]), [rows, sort]);
  return <div className="table-scroll compact-table"><table><thead><tr><SortHeader label="Input statistic" column="label" state={sort} onChange={setSort} /><SortHeader label="Reads · n" column="count" state={sort} onChange={setSort} /><SortHeader label="% of FASTQ records" column="percent" state={sort} onChange={setSort} /></tr></thead><tbody>{displayed.map((row) => <tr key={row.key}><td>{row.label}{row.note && <small className="cell-note">{row.note}</small>}</td><td>{row.count.toLocaleString()}</td><td>{pct(row.percent)}</td></tr>)}</tbody></table></div>;
}

function ParameterSettingsTable({ bundle }: { bundle: ResultBundle }) {
  const [sort, setSort] = useState<SortState<string>>({ key: "scope", direction: "asc" });
  const rows = useMemo(() => sorted(parameterSettings(bundle), sort, (row, key) => row[key as keyof ParameterSettingRow]), [bundle, sort]);
  return <article className="statistics-card run-wide-stats"><header><h3>Parameters and effective sample settings</h3><p>Every stored run parameter is shown, followed by each sample’s inputs and effective override values.</p></header><div className="table-scroll parameter-table"><table><thead><tr><SortHeader label="Scope" column="scope" state={sort} onChange={setSort} /><SortHeader label="Sample" column="sample" state={sort} onChange={setSort} /><SortHeader label="Parameter" column="parameter" state={sort} onChange={setSort} /><SortHeader label="Value used" column="value" state={sort} onChange={setSort} /></tr></thead><tbody>{rows.map((row, index) => <tr key={`${row.scope}-${row.sample}-${row.parameter}-${index}`}><td>{row.scope}</td><td>{row.sample}</td><td><code>{row.parameter}</code></td><td><code title={row.value}>{row.value}</code></td></tr>)}</tbody></table></div></article>;
}

const OPTIONAL_STAGE_LABELS: Record<OptionalStageName, string> = {
  contamination: "Contamination checks", postprocessing: "Alignment + downstream filtering",
  collapse: "Haplotype collapse + functional filtering", tree: "Collapsed phylogenies",
};

function OptionalStagePanel({ bundle, progress, onCompute, onSkip }: { bundle: ResultBundle; progress?: DeferredAnalysisProgress;
  onCompute(stage: OptionalStageName): void; onSkip(): void }) {
  return <section className="optional-stage-panel"><header><div><span className="section-kicker">Resumable analysis</span><h2>Optional stages after consensus</h2><p>Deferred and skipped work is recorded explicitly. Contamination is an optional side gate: downstream analysis can continue unfiltered. Collapsed-variant functional filtering runs with collapse, which requires post-processing; the default tree then uses those collapsed variants.</p><MethodLink topic="optional-stages" /></div></header>
    <div className="optional-stage-grid">{OPTIONAL_STAGE_ORDER.map((stage) => {
      const status = stageStatus(bundle, stage), prerequisites: OptionalStageName[] = stage === "collapse" ? ["postprocessing"] : stage === "tree" ? ["postprocessing", "collapse"] : [];
      const missing = prerequisites.filter((candidate) => !stageCompleted(bundle, candidate));
      const needsReapply = stage === "postprocessing" && status.state === "completed" && stageCompleted(bundle, "contamination")
        && bundle.postprocessingContaminationMode === "bypassed" && bundle.config.parameters.contaminationFilter;
      const busy = Boolean(progress), current = progress?.stage === stage;
      return <article className={`optional-stage ${status.state}${needsReapply ? " needs-reapply" : ""}`} key={stage}><span>{OPTIONAL_STAGE_LABELS[stage]}</span><strong>{needsReapply ? "computed · unfiltered" : status.state}</strong><p>{status.detail}</p>{(status.state !== "completed" || needsReapply) && <button type="button" disabled={busy} onClick={() => onCompute(stage)}>{needsReapply ? "Recompute and apply contamination" : `Compute ${OPTIONAL_STAGE_LABELS[stage].toLowerCase()}`}</button>}{missing.length > 0 && <small>Also computes required input: {missing.map((item) => OPTIONAL_STAGE_LABELS[item]).join(" → ")}</small>}{current && <div className="on-demand-progress"><progress max="100" value={Math.round(progress.fraction * 100)} /><span>{progress.detail}</span><button type="button" className="danger" onClick={onSkip}>Skip this step</button></div>}</article>;
    })}</div>
  </section>;
}

function AllSampleOverview({ bundle, onOpenSample }: { bundle: ResultBundle; onOpenSample(sample: string): void }) {
  const [sort, setSort] = useState<SortState<string>>({ key: "sample", direction: "asc" });
  const rows = useMemo(() => sorted(sampleOverviewStats(bundle), sort, (row, key) => row[key as keyof SampleOverviewStat]), [bundle, sort]);
  const percentCell = (value: number) => <span className={value > 0 ? "stat-reject" : ""}>{pct(value)}</span>;
  const contaminationComplete = stageCompleted(bundle, "contamination"), postprocessingComplete = stageCompleted(bundle, "postprocessing"), collapseComplete = stageCompleted(bundle, "collapse");
  return <section className="result-section all-sample-overview"><div className="table-heading"><div><span className="section-kicker">Across-sample overview</span><h2>Every configured sample</h2><p>Subsampling percentages use demultiplexed reads. UMI-filter percentages use selected reads represented by stored family calls; consensus-filter percentages use reads represented by consensus families.</p><MethodLink topic="results" /></div><span>{rows.length.toLocaleString()} samples</span></div><div className="table-scroll overview-table"><table><thead><tr>
    <SortHeader label="Sample" column="sample" state={sort} onChange={setSort} /><SortHeader label="Donor ID" column="donorId" state={sort} onChange={setSort} /><SortHeader label="Demux reads" column="demultiplexedReads" state={sort} onChange={setSort} /><SortHeader label="Selected reads" column="selectedReads" state={sort} onChange={setSort} /><SortHeader label="Subsampled" column="downsampledPercent" state={sort} onChange={setSort} /><SortHeader label="Observed UMI families" column="observedFamilies" state={sort} onChange={setSort} /><SortHeader label="BPB read rejects" column="bpbReadPercent" state={sort} onChange={setSort} /><SortHeader label="UMI-length read rejects" column="umiLengthReadPercent" state={sort} onChange={setSort} /><SortHeader label="Family-size read rejects" column="familySizeReadPercent" state={sort} onChange={setSort} /><SortHeader label="LDA read rejects" column="ldaReadPercent" state={sort} onChange={setSort} /><SortHeader label="Heteroduplex read rejects" column="heteroduplexReadPercent" state={sort} onChange={setSort} /><SortHeader label="Consensus families" column="consensusFamilies" state={sort} onChange={setSort} /><SortHeader label="Artefact read rejects" column="artefactReadPercent" state={sort} onChange={setSort} /><SortHeader label="Agreement read rejects" column="agreementReadPercent" state={sort} onChange={setSort} /><SortHeader label="Contam read rejects" column="contaminationReadPercent" state={sort} onChange={setSort} /><SortHeader label="Panel read rejects" column="panelReadPercent" state={sort} onChange={setSort} /><SortHeader label="Functional represented-family rejects" column="functionalReadPercent" state={sort} onChange={setSort} /><SortHeader label="Retained families" column="retainedFamilies" state={sort} onChange={setSort} /><SortHeader label="Functional variants evaluated" column="functionalEvaluatedFamilies" state={sort} onChange={setSort} /><SortHeader label="Functional variants passed" column="functionalPassedFamilies" state={sort} onChange={setSort} /><SortHeader label="Collapsed haplotypes" column="collapsedHaplotypes" state={sort} onChange={setSort} />
  </tr></thead><tbody>{rows.map((row) => <tr key={row.sample}><td><button type="button" className="sample-link" onClick={() => onOpenSample(row.sample)}>{row.sample}</button></td><td>{row.donorId}</td><td>{row.demultiplexedReads.toLocaleString()}</td><td>{row.selectedReads.toLocaleString()}</td><td>{row.downsampledReads.toLocaleString()} <small>({pct(row.downsampledPercent)})</small></td><td>{row.observedFamilies.toLocaleString()}</td><td>{percentCell(row.bpbReadPercent)}</td><td>{percentCell(row.umiLengthReadPercent)}</td><td>{percentCell(row.familySizeReadPercent)}</td><td>{percentCell(row.ldaReadPercent)}</td><td>{percentCell(row.heteroduplexReadPercent)}</td><td>{row.consensusFamilies.toLocaleString()}</td><td>{postprocessingComplete ? percentCell(row.artefactReadPercent) : "not computed"}</td><td>{postprocessingComplete ? percentCell(row.agreementReadPercent) : "not computed"}</td><td>{contaminationComplete ? percentCell(row.contaminationReadPercent) : "not computed"}</td><td>{postprocessingComplete ? percentCell(row.panelReadPercent) : "not computed"}</td><td>{!row.functionalConfigured ? "—" : collapseComplete ? percentCell(row.functionalReadPercent) : "not computed"}</td><td>{postprocessingComplete ? row.retainedFamilies.toLocaleString() : "not computed"}</td><td>{!row.functionalConfigured ? "—" : collapseComplete ? row.functionalEvaluatedFamilies.toLocaleString() : "not computed"}</td><td>{!row.functionalConfigured ? "—" : collapseComplete ? row.functionalPassedFamilies.toLocaleString() : "not computed"}</td><td>{collapseComplete ? row.collapsedHaplotypes.toLocaleString() : "not computed"}</td></tr>)}</tbody></table></div>
    <article className="statistics-card run-wide-stats"><header><h3>Run-wide input filtering</h3><p>These filters occur before sample assignment, so they cannot be attributed scientifically to individual samples. Percentages use all FASTQ records as the denominator.</p></header><CountStatsTable rows={inputFilterStats(bundle)} /></article>
    <ParameterSettingsTable bundle={bundle} />
    <TimingSummary bundle={bundle} />
  </section>;
}

function changeDetails(changes: AlignmentChangeSummary | undefined): string[] {
  if (!changes) return [];
  const details = [`${changes.rowsBefore} → ${changes.rowsAfter} rows`, `${changes.columnsBefore} → ${changes.columnsAfter} columns`];
  if (changes.rowOrderChanged) details.push(`Row order changed: [${changes.rowOrderBefore.join(", ")}] → [${changes.rowOrderAfter.join(", ")}]`);
  if (changes.removedRows.length) details.push(`Deleted rows: ${changes.removedRows.join(", ")}`);
  if (changes.addedRows.length) details.push(`Added/renamed rows: ${changes.addedRows.join(", ")}`);
  if (changes.changedRows.length) details.push(`Changed rows: ${changes.changedRows.join(", ")}`);
  for (const row of changes.rowChanges ?? []) {
    const operations = [`${row.substitutedNucleotides} substitutions`, `${row.insertedNucleotides} inserted bases`, `${row.removedNucleotides} deleted bases`];
    if (row.gapPlacementChanged) operations.push("gap placement changed with identical ungapped sequence");
    details.push(`${row.name}: ${operations.join(", ")}`);
  }
  if (changes.substitutedNucleotides) details.push(`${changes.substitutedNucleotides} substituted bases`);
  if (changes.insertedNucleotides) details.push(`${changes.insertedNucleotides} inserted bases`);
  if (changes.removedNucleotides) details.push(`${changes.removedNucleotides} deleted bases`);
  if (!details.slice(2).length) details.push("Gap placement and/or translation settings only; ungapped nucleotide content is unchanged");
  return details;
}

function appendAlignmentAudit(bundle: ResultBundle, entry: AlignmentAuditEntry): ResultBundle {
  return { ...bundle, alignmentEditHistory: [...(bundle.alignmentEditHistory ?? []), entry] };
}

function AlignmentAuditTable({ bundle, sample }: { bundle: ResultBundle; sample: string }) {
  const [sort, setSort] = useState<SortState<string>>({ key: "timestamp", direction: "desc" });
  const rows = useMemo(() => {
    const stored = (bundle.alignmentEditHistory ?? []).filter((entry) => entry.alignmentKey.startsWith(`${sample}/`));
    if (stored.length) return stored;
    return Object.entries(bundle.alignmentEdits ?? {}).filter(([key]) => key.startsWith(`${sample}/`)).map(([alignmentKey, edit]) => ({
      alignmentKey, action: "alignment-edit" as const, timestamp: edit.savedUtc, source: edit.source,
      details: changeDetails(edit.changes).length ? changeDetails(edit.changes) : edit.warnings ?? ["Legacy edit; detailed change counts were not stored."],
      beforeFingerprint: edit.baselineFingerprint, afterFingerprint: edit.editedFingerprint,
    }));
  }, [bundle, sample]);
  const displayed = useMemo(() => sorted(rows, sort, (row, key) => key === "details" ? row.details.join("; ") : row[key as keyof AlignmentAuditEntry]), [rows, sort]);
  if (!displayed.length) return <p className="no-audit-edits">No interactive alignment edits or tree recalculations have been recorded for this sample.</p>;
  return <div className="table-scroll compact-table audit-table"><table><thead><tr><SortHeader label="Time" column="timestamp" state={sort} onChange={setSort} /><SortHeader label="Alignment" column="alignmentKey" state={sort} onChange={setSort} /><SortHeader label="Action" column="action" state={sort} onChange={setSort} /><SortHeader label="Source" column="source" state={sort} onChange={setSort} /><SortHeader label="Exact recorded changes" column="details" state={sort} onChange={setSort} /></tr></thead><tbody>{displayed.map((row, index) => <tr key={`${row.timestamp}-${row.alignmentKey}-${index}`}><td>{row.timestamp}</td><td><code>{row.alignmentKey}</code></td><td>{row.action.replaceAll("-", " ")}</td><td>{row.source}</td><td>{row.details.map((detail) => <span className="audit-detail" key={detail}>{detail}</span>)}</td></tr>)}</tbody></table></div>;
}

function ThresholdAuditTrail({ bundle, sample }: { bundle: ResultBundle; sample: string }) {
  const records = bundle.thresholdSelections ?? [];
  const samplePrefixes = bundle.config.samples.map((row) => `${row.name}.`);
  if (!records.length) return <p className="no-audit-edits">Interactive filtering was not enabled for this run.</p>;
  return <ol className="threshold-audit">{records.map((record) => <li key={record.id}><strong>{record.phase === "umi" ? "UMI-family decision" : "Consensus-filter decision"}</strong><time>{record.acceptedUtc}</time><ul>{record.changes.filter((change) => !samplePrefixes.some((prefix) => change.startsWith(prefix)) || change.startsWith(`${sample}.`)).map((change) => <li key={change}>{change}</li>)}</ul></li>)}</ol>;
}

function editWarnings(inspection: ReturnType<typeof validateCorrectedAlignment>): string[] {
  const warnings: string[] = [];
  if (inspection.removedRows.length) warnings.push(`${inspection.removedRows.length} sequence rows deleted`);
  if (inspection.addedRows.length) warnings.push(`${inspection.addedRows.length} sequence rows added or renamed`);
  if (inspection.substitutedNucleotides) warnings.push(`${inspection.substitutedNucleotides} base positions substituted`);
  if (inspection.insertedNucleotides) warnings.push(`${inspection.insertedNucleotides} nucleotide characters inserted`);
  if (inspection.removedNucleotides) warnings.push(`${inspection.removedNucleotides} nucleotide characters deleted`);
  return warnings;
}

function referenceSequence(bundle: ResultBundle, sample: string, variant: AlignmentVariant = "collapsed"): string | undefined {
  const fasta = bundle.referenceAlignments?.[alignmentKey(sample, variant)]
    ?? (variant === "uncollapsed" ? bundle.referenceAlignments?.[alignmentKey(sample)] : undefined);
  if (!fasta) return undefined;
  try { return parseFasta(fasta)[0]?.sequence; } catch { return undefined; }
}

function collapsedMetadata(groups: CollapseGroup[] | undefined): Record<string, LeafMetadata> {
  return Object.fromEntries((groups ?? []).map((group) => [group.representativeId, { familyCount: group.familyCount }]));
}

function uncollapsedMetadata(records: PostprocRecord[]): Record<string, LeafMetadata> {
  return Object.fromEntries(records.filter((record) => record.alignedNt).map((record) => [record.id, { familyCount: 1, minimumAgreement: record.minimumAgreement }]));
}

const DONOR_COLORS = ["#08796f", "#c5534f", "#5e55a4", "#c77b20", "#3c78a8", "#8b5a82", "#557b3d", "#a84b70"];

function DonorView({ bundle, donorId }: { bundle: ResultBundle; donorId: string }) {
  const [variant, setVariant] = useState<"collapsed" | "functional">("collapsed");
  const [analysis, setAnalysis] = useState<DonorPhylogeny>();
  const [busy, setBusy] = useState(false), [error, setError] = useState("");
  const [alphabet, setAlphabet] = useState<"nt" | "aa">("nt");
  const samples = useMemo(() => bundle.config.samples.filter((sample) => sample.donorId === donorId), [bundle, donorId]);
  const summaries = useMemo(() => samples.map((configured) => {
    const summary = bundle.summaries.find((row) => row.sample === configured.name);
    return { sample: configured.name, consensuses: summary?.consensusSequences ?? 0,
      retained: summary?.postprocPassed ?? 0, collapsed: bundle.collapseGroups?.[configured.name]?.length ?? 0,
      functional: summary?.functionalPassed ?? 0 };
  }), [bundle, samples]);
  useEffect(() => { setAnalysis(undefined); setError(""); setVariant("collapsed"); }, [donorId]);

  async function build() {
    setBusy(true); setError("");
    try {
      const rows: Array<{ name: string; sequence: string; sample: string; familyCount: number }> = [];
      for (const configured of samples) {
        const active = effectiveAlignment(bundle, configured.name, variant).fasta; if (!active) continue;
        const familyCounts = new Map((bundle.collapseGroups?.[configured.name] ?? []).map((group) => [group.representativeId, group.familyCount]));
        parseFasta(active).forEach((record, index) => rows.push({
          name: `${configured.name.replace(/[^A-Za-z0-9_.+-]/g, "_")}__${index + 1}__${record.name}`,
          sequence: record.sequence.replaceAll("-", "").toUpperCase().replaceAll("U", "T"), sample: configured.name,
          familyCount: familyCounts.get(record.name) ?? 1,
        }));
      }
      if (!rows.length) throw new Error(`No ${variant} sequences are available for donor ${donorId}. Compute the required downstream stages first.`);
      const aligned = rows.length > 1 ? await runScalableMsa(rows.map((row) => row.sequence), runAlivibeMsa, undefined, 3, "nucleotide") : [rows[0].sequence];
      const fasta = rows.map((row, index) => `>${row.name}\n${aligned[index]}\n`).join(""), newick = await runFastTree(fasta);
      const legend = samples.map((sample, index) => ({ label: sample.name, color: DONOR_COLORS[index % DONOR_COLORS.length] }));
      const colors = new Map(legend.map((row) => [row.label, row.color]));
      setAnalysis({ donorId, variant, fasta, newick, legend,
        leafMetadata: Object.fromEntries(rows.map((row) => [row.name, { familyCount: row.familyCount,
          color: colors.get(row.sample), category: row.sample }])) });
    } catch (cause) { setError(cause instanceof Error ? cause.message : String(cause)); }
    finally { setBusy(false); }
  }

  return <section className="result-section donor-view"><div className="table-heading"><div><span className="section-kicker">Donor-level view</span><h2>{donorId}</h2><p>Combined sample sequences retain their sample identity. Tip color records the source sample; collapsed bubble area remains exactly linear in represented UMI-family count.</p><MethodLink topic="contamination" /></div><span>{samples.length} sample{samples.length === 1 ? "" : "s"}</span></div>
    <div className="donor-sample-grid">{summaries.map((row) => <article key={row.sample}><strong>{row.sample}</strong><dl><div><dt>Consensus families</dt><dd>{row.consensuses.toLocaleString()}</dd></div><div><dt>Postproc retained</dt><dd>{row.retained.toLocaleString()}</dd></div><div><dt>Collapsed haplotypes</dt><dd>{row.collapsed.toLocaleString()}</dd></div><div><dt>Functional variants</dt><dd>{row.functional.toLocaleString()}</dd></div></dl></article>)}</div>
    <section className="donor-workbench"><header><div><span className="section-kicker">Joint alignment + phylogeny</span><h3>{variant === "collapsed" ? "Post-collapse haplotypes" : "Functional collapsed variants"}</h3><p>A new joint MSA is calculated across all samples in this donor before bundled double-precision FastTree inference.</p></div><div className="donor-variant"><button type="button" className={variant === "collapsed" ? "active" : ""} onClick={() => { setVariant("collapsed"); setAnalysis(undefined); }}>Collapsed</button><button type="button" className={variant === "functional" ? "active" : ""} onClick={() => { setVariant("functional"); setAnalysis(undefined); }}>Functional</button><button type="button" className="primary" disabled={busy} onClick={() => void build()}>{busy ? "Aligning + inferring…" : analysis?.variant === variant ? "Rebuild" : "Build tree + alignment"}</button></div></header>
      {error && <div className="error-box" role="alert">{error}</div>}
      {analysis?.donorId === donorId && analysis.variant === variant && <AlignmentTreeViewer fasta={analysis.fasta} newick={analysis.newick} alphabet={alphabet} onAlphabetChange={setAlphabet} leafMetadata={analysis.leafMetadata} tipLegend={analysis.legend} collapsed variantLabel={`Donor ${donorId} · ${variant}`} name={`${safeDatasetName(bundle.config.dataset)}-donor-${safeDatasetName(donorId)}-${variant}`} />}
    </section>
  </section>;
}

export function ResultsExplorer({ bundle, onSaveResults, onBundleChange }: { bundle: ResultBundle; onSaveResults(): void; onBundleChange(bundle: ResultBundle): void }) {
  const [tab, setTab] = useState<Tab>("overview"), [sample, setSample] = useState(bundle.summaries[0]?.sample ?? bundle.config.samples[0]?.name ?? ""), [allSamples, setAllSamples] = useState(true);
  const [donorView, setDonorView] = useState("");
  const [query, setQuery] = useState(""), [familyQuery, setFamilyQuery] = useState(""), [alphabet, setAlphabet] = useState<"nt" | "aa">("nt"), [showUncollapsed, setShowUncollapsed] = useState(false), [showFunctional, setShowFunctional] = useState(false);
  const [familyPage, setFamilyPage] = useState(0), [sequencePage, setSequencePage] = useState(0), [contaminationPage, setContaminationPage] = useState(0);
  const [familySort, setFamilySort] = useState<SortState<string>>({ key: "familySize", direction: "desc" });
  const [sequenceSort, setSequenceSort] = useState<SortState<string>>({ key: "minimumAgreement", direction: "asc" });
  const [contaminationSort, setContaminationSort] = useState<SortState<string>>({ key: "nearestNonselfDistance", direction: "asc" });
  const [alignmentStatus, setAlignmentStatus] = useState(""), [alignmentError, setAlignmentError] = useState(""), [alignmentBusy, setAlignmentBusy] = useState("");
  const [contaminationPhylogeny, setContaminationPhylogeny] = useState<ContaminationPhylogeny>(), [contaminationPhylogenyBusy, setContaminationPhylogenyBusy] = useState(false), [contaminationPhylogenyError, setContaminationPhylogenyError] = useState("");
  const [deferredProgress, setDeferredProgress] = useState<DeferredAnalysisProgress>(), [deferredError, setDeferredError] = useState("");
  const bundleRef = useRef(bundle), sampleRef = useRef(sample), changeRef = useRef(onBundleChange), alivibeRef = useRef<AlivibeSession | null>(null);
  const deferredControllerRef = useRef<AbortController | undefined>(undefined), deferredStageRef = useRef<OptionalStageName | undefined>(undefined);
  bundleRef.current = bundle; sampleRef.current = sample; changeRef.current = onBundleChange;

  const summary = bundle.summaries.find((row) => row.sample === sample), sampleConfig = bundle.config.samples.find((row) => row.name === sample);
  const families = useMemo(() => bundle.umiFamilies.filter((row) => row.sample === sample), [bundle, sample]);
  const consensuses = useMemo(() => bundle.consensuses.filter((row) => row.sample === sample), [bundle, sample]);
  const currentOverview = useMemo(() => sampleOverviewStats(bundle).find((row) => row.sample === sample), [bundle, sample]);
  const porpidStats = useMemo(() => porpidCallStats(bundle, sample), [bundle, sample]);
  const postprocStats = useMemo(() => postprocFilterStats(bundle, sample), [bundle, sample]);
  const functionalStats = useMemo(() => functionalFilterStats(bundle, sample), [bundle, sample]);
  const filteredFamilies = useMemo(() => sorted(families.filter((row) => !familyQuery || row.umi.includes(familyQuery.toUpperCase()) || row.mostLikelyParent.includes(familyQuery.toUpperCase()) || row.disposition.toLowerCase().includes(familyQuery.toLowerCase())), familySort,
    (row: UmiFamily, key) => key === "classification" ? row.disposition : row[key as keyof UmiFamily]), [families, familyQuery, familySort]);
  const sampleRecords = useMemo(() => bundle.records.filter((row) => row.sample === sample), [bundle, sample]);
  const records = useMemo(() => sorted(sampleRecords.filter((row) => !query || row.id.toLowerCase().includes(query.toLowerCase()) || row.umi.includes(query.toUpperCase())), sequenceSort,
    (row: PostprocRecord, key) => key === "filters" ? Number(row.artefactPass && row.agreementPass && row.contaminationPass && row.panelPass) : key === "apobec" ? row.apobec?.posteriorGaInflated : key === "reason" ? row.rejectionReasons.join(";") : row[key as keyof PostprocRecord]), [sampleRecords, query, sequenceSort]);
  const contaminationRows = useMemo(() => sorted(deduplicateContaminationCalls(bundle.contamination).filter((row) => row.sample === sample), contaminationSort,
    (row: ContaminationCall, key) => key === "decision" ? row.discarded ? 2 : row.suspectOnly ? 1 : 0 : row[key as keyof ContaminationCall]), [bundle, sample, contaminationSort]);
  const collapsed = useMemo(() => effectiveAlignment(bundle, sample, "collapsed"), [bundle, sample]);
  const uncollapsed = useMemo(() => effectiveAlignment(bundle, sample, "uncollapsed"), [bundle, sample]);
  const functional = useMemo(() => effectiveAlignment(bundle, sample, "functional"), [bundle, sample]);
  const collapsedTree = collapsed.edit?.treeNewick ?? bundle.trees[collapsed.key], uncollapsedTree = uncollapsed.edit?.treeNewick ?? bundle.trees[uncollapsed.key];
  const functionalTree = functional.edit?.treeNewick ?? bundle.trees[functional.key];
  const refSequence = useMemo(() => referenceSequence(bundle, sample), [bundle, sample]);
  const functionalRefSequence = useMemo(() => referenceSequence(bundle, sample, "functional"), [bundle, sample]);
  const collapsedTips = useMemo(() => collapsedMetadata(bundle.collapseGroups?.[sample]), [bundle, sample]);
  const uncollapsedTips = useMemo(() => uncollapsedMetadata(sampleRecords), [sampleRecords]);
  const functionalTips = useMemo(() => collapsedMetadata(bundle.collapseGroups?.[sample]?.filter((group) => group.functionalPass === true)), [bundle, sample]);
  const edited = Boolean(collapsed.edit || uncollapsed.edit || functional.edit);
  const contaminationDone = stageCompleted(bundle, "contamination"), postprocessingDone = stageCompleted(bundle, "postprocessing"), collapseDone = stageCompleted(bundle, "collapse");
  const donorIds = useMemo(() => [...new Set(bundle.config.samples.map((row) => row.donorId).filter((value): value is string => Boolean(value)))].sort(), [bundle]);
  const sampleMode = !allSamples && !donorView;

  useEffect(() => { setFamilyPage(0); setSequencePage(0); setContaminationPage(0); setAlignmentStatus(""); setAlignmentError(""); setShowUncollapsed(false); setShowFunctional(false); setContaminationPhylogeny(undefined); setContaminationPhylogenyError(""); }, [sample]);
  useEffect(() => setFamilyPage(0), [familyQuery, familySort]); useEffect(() => setSequencePage(0), [query, sequenceSort]); useEffect(() => setContaminationPage(0), [contaminationSort]);
  useEffect(() => () => { deferredControllerRef.current?.abort(); const session = alivibeRef.current; if (session?.timer) window.clearInterval(session.timer); if (session && !session.popup.closed) session.popup.close(); }, []);
  const page = <T,>(rows: T[], index: number) => rows.slice(index * PAGE_SIZE, (index + 1) * PAGE_SIZE);

  async function computeOptionalStage(target: OptionalStageName) {
    if (deferredControllerRef.current) return;
    const currentBundle = bundleRef.current;
    const reapplyingContamination = target === "postprocessing" && stageCompleted(currentBundle, "contamination")
      && currentBundle.postprocessingContaminationMode === "bypassed" && currentBundle.config.parameters.contaminationFilter;
    if (reapplyingContamination && !window.confirm("Recomputing downstream filtering will apply the newly available contamination decisions, replace downstream alignments, invalidate collapsed trees, and detach any active alignment edits. The original consensus calls remain unchanged. Continue?")) return;
    const controller = new AbortController(); deferredControllerRef.current = controller; deferredStageRef.current = target;
    setDeferredError(""); setDeferredProgress({ stage: target, fraction: 0, detail: `Preparing ${OPTIONAL_STAGE_LABELS[target].toLowerCase()}…` });
    try {
      const next = await computeThrough(bundleRef.current, target, { signal: controller.signal,
        onProgress: (state) => { deferredStageRef.current = state.stage; setDeferredProgress(state); },
        onCheckpoint: (checkpointBundle) => { bundleRef.current = checkpointBundle; changeRef.current(checkpointBundle); },
      });
      bundleRef.current = next; changeRef.current(next);
    } catch (cause) {
      if (controller.signal.aborted) {
        const skipped = markOptionalStageSkipped(bundleRef.current, deferredStageRef.current ?? target);
        bundleRef.current = skipped; changeRef.current(skipped);
      } else setDeferredError(cause instanceof Error ? cause.message : String(cause));
    } finally {
      if (deferredControllerRef.current === controller) deferredControllerRef.current = undefined;
      deferredStageRef.current = undefined; setDeferredProgress(undefined);
    }
  }

  function skipDeferredStage() { deferredControllerRef.current?.abort(); }

  async function installCorrection(variant: AlignmentVariant, targetSample: string, baseline: string, corrected: string, frameOffset: AlignmentFrameOffset, source: string) {
    setAlignmentBusy(`${targetSample}/${variant}`); setAlignmentError("");
    try {
      const inspected = validateCorrectedAlignment(baseline, corrected), confirmationWarnings = editWarnings(inspected);
      if (confirmationWarnings.length && !window.confirm(`This alignment edit changes biological content:\n\n• ${confirmationWarnings.join("\n• ")}\n\nThe pipeline-generated alignment will remain preserved. Save this separate edit?`)) {
        setAlignmentStatus("Alignment edit cancelled; the stored alignment was not changed."); return;
      }
      const currentBundle = bundleRef.current, key = alignmentKey(targetSample, variant), original = currentBundle.alignments[key] ?? (variant === "uncollapsed" ? currentBundle.alignments[alignmentKey(targetSample)] : undefined);
      if (!original) throw new Error("The original nucleotide alignment is unavailable.");
      const current = effectiveAlignment(currentBundle, targetSample, variant), existingTree = current.edit?.treeNewick ?? currentBundle.trees[key];
      const originalInspection = validateCorrectedAlignment(original, inspected.fasta), warnings = editWarnings(originalInspection);
      const changes = summarizeAlignmentChanges(original, inspected.fasta), incrementalChanges = summarizeAlignmentChanges(baseline, inspected.fasta), timestamp = new Date().toISOString();
      let next: ResultBundle = { ...currentBundle, alignmentEdits: { ...currentBundle.alignmentEdits, [key]: {
        fasta: inspected.fasta, frameOffset, baselineFingerprint: inspectAlignment(original, 1).fingerprint, editedFingerprint: inspected.fingerprint,
        source, savedUtc: timestamp, treeNewick: existingTree,
        treeFingerprint: current.edit?.treeFingerprint ?? (existingTree ? inspectAlignment(current.fasta ?? original, 1).fingerprint : undefined),
        treeStale: Boolean(existingTree) && inspected.fingerprint !== (current.edit?.treeFingerprint ?? inspectAlignment(current.fasta ?? original, 1).fingerprint), warnings, changes,
      } }, log: [...currentBundle.log, `${timestamp} alignment edit: ${key}; ${changeDetails(incrementalChanges).join("; ")}; source=${source}`] };
      next = appendAlignmentAudit(next, { alignmentKey: key, action: "alignment-edit", timestamp, source,
        details: changeDetails(incrementalChanges), beforeFingerprint: inspectAlignment(baseline, 1).fingerprint, afterFingerprint: inspected.fingerprint });
      bundleRef.current = next; changeRef.current(next);
      setAlignmentStatus(`Saved a separate ${variant} alignment edit: ${inspected.rows.toLocaleString()} rows × ${inspected.columns.toLocaleString()} columns.${warnings.length ? ` Warning: ${warnings.join("; ")}.` : " Gap placement changed; nucleotide content is unchanged."} Recalculate the tree when ready.`);
    } catch (cause) { setAlignmentError(cause instanceof Error ? cause.message : String(cause)); }
    finally { setAlignmentBusy(""); }
  }

  function setFrameOffset(variant: AlignmentVariant, frameOffset: AlignmentFrameOffset) {
    const currentBundle = bundleRef.current, targetSample = sampleRef.current, current = effectiveAlignment(currentBundle, targetSample, variant), key = current.key;
    if (!current.fasta) return;
    if (current.frameOffset === frameOffset) return;
    const fingerprint = inspectAlignment(current.fasta, 1).fingerprint, original = currentBundle.alignments[key] ?? current.fasta, timestamp = new Date().toISOString();
    let next: ResultBundle = { ...currentBundle, alignmentEdits: { ...currentBundle.alignmentEdits, [key]: {
      fasta: current.fasta, frameOffset, baselineFingerprint: inspectAlignment(original, 1).fingerprint, editedFingerprint: fingerprint,
      source: current.edit?.source ?? "Translation frame selection", savedUtc: timestamp, treeNewick: current.edit?.treeNewick ?? currentBundle.trees[key],
      treeFingerprint: current.edit?.treeFingerprint ?? (currentBundle.trees[key] ? fingerprint : undefined), treeStale: current.edit?.treeStale ?? false, warnings: current.edit?.warnings,
      changes: current.edit?.changes ?? summarizeAlignmentChanges(original, current.fasta),
    } }, log: [...currentBundle.log, `${timestamp} translation frame: ${key}; nucleotide column ${current.frameOffset + 1} -> ${frameOffset + 1}`] };
    next = appendAlignmentAudit(next, { alignmentKey: key, action: "frame-change", timestamp, source: "Translation frame selection",
      details: [`Protein translation start changed from nucleotide column ${current.frameOffset + 1} to ${frameOffset + 1}`], beforeFingerprint: fingerprint, afterFingerprint: fingerprint });
    bundleRef.current = next; changeRef.current(next); setAlignmentStatus(`Protein translation now starts at nucleotide column ${frameOffset + 1}; this choice is stored without changing the original alignment.`);
  }

  function resetCorrection(variant: AlignmentVariant) {
    const currentBundle = bundleRef.current, key = alignmentKey(sampleRef.current, variant), previous = currentBundle.alignmentEdits?.[key], edits = { ...currentBundle.alignmentEdits }, timestamp = new Date().toISOString();
    delete edits[key]; let next: ResultBundle = { ...currentBundle, alignmentEdits: Object.keys(edits).length ? edits : undefined,
      log: [...currentBundle.log, `${timestamp} alignment reset: ${key}; restored immutable pipeline alignment`] };
    next = appendAlignmentAudit(next, { alignmentKey: key, action: "edit-reset", timestamp, source: "Results workbench",
      details: ["Discarded the active edited copy and restored the pipeline alignment", ...changeDetails(previous?.changes)],
      beforeFingerprint: previous?.editedFingerprint, afterFingerprint: previous?.baselineFingerprint });
    bundleRef.current = next; changeRef.current(next); setAlignmentError(""); setAlignmentStatus(`Restored the pipeline-generated ${variant} alignment and its stored tree.`);
  }

  async function inferTree(variant: AlignmentVariant) {
    const currentBundle = bundleRef.current, targetSample = sampleRef.current, current = effectiveAlignment(currentBundle, targetSample, variant), key = current.key;
    if (!current.fasta) return;
    setAlignmentBusy(`${targetSample}/${variant}`); setAlignmentError(""); setAlignmentStatus(`Inferring the ${variant} tree with bundled double-precision FastTree…`);
    try {
      const treeNewick = await runFastTree(current.fasta), fingerprint = inspectAlignment(current.fasta, 1).fingerprint;
      let next: ResultBundle;
      const timestamp = new Date().toISOString();
      if (current.edit) next = { ...currentBundle, alignmentEdits: { ...currentBundle.alignmentEdits, [key]: { ...current.edit, treeNewick, treeFingerprint: fingerprint, treeStale: false, savedUtc: timestamp } }, log: [...currentBundle.log, `${timestamp} interactive FastTree: ${key}; edited alignment fingerprint ${fingerprint}`] };
      else next = { ...currentBundle, trees: { ...currentBundle.trees, [key]: treeNewick }, log: [...currentBundle.log, `${timestamp} interactive FastTree: ${key}; alignment fingerprint ${fingerprint}`] };
      next = appendAlignmentAudit(next, { alignmentKey: key, action: "tree-recalculation", timestamp, source: "Bundled double-precision FastTree",
        details: [`Recalculated ${variant} topology and branch lengths from the active alignment`, `Alignment fingerprint: ${fingerprint}`], beforeFingerprint: fingerprint, afterFingerprint: fingerprint });
      const label = variant === "collapsed" ? "Collapsed" : variant === "uncollapsed" ? "Uncollapsed" : "Functional-sequence";
      bundleRef.current = next; changeRef.current(next); setAlignmentStatus(`${label} tree inference complete. Mutation mapping has been refreshed for this alignment.`);
    } catch (cause) { setAlignmentError(cause instanceof Error ? cause.message : String(cause)); }
    finally { setAlignmentBusy(""); }
  }

  function createMsaJob(sequences: string[], scoringMode: "nucleotide" | "amino-acid"): AlivibeMsaJob {
    const controller = new AbortController(); return { result: runAlivibeMsa(sequences, controller.signal, 3, scoringMode), cancel: () => controller.abort() };
  }

  async function buildContaminationPhylogeny() {
    const currentBundle = bundleRef.current, targetSample = sampleRef.current;
    setContaminationPhylogenyBusy(true); setContaminationPhylogenyError("");
    try {
      const references = currentBundle.contaminationReferences ?? [];
      if (!references.length) throw new Error("This results file does not contain the contamination-panel reference sequences. Re-run with webPORPID 0.3.5 or later to build this view.");
      const donor = currentBundle.consensuses.filter((record) => record.sample === targetSample);
      if (!donor.length) throw new Error("This sample has no consensus sequences to place beside the contamination references.");
      const decisions = new Map(deduplicateContaminationCalls(currentBundle.contamination).filter((call) => call.sample === targetSample).map((call) => [call.sequenceId, call]));
      const rows = [
        ...references.map((record, index) => ({ name: `reference_${index + 1}__${record.name}`, sequence: record.sequence.replaceAll("-", "").toUpperCase().replaceAll("U", "T"), color: CONTAMINATION_TIP_LEGEND[0].color, category: CONTAMINATION_TIP_LEGEND[0].label })),
        ...donor.map((record) => { const contaminant = Boolean(decisions.get(record.id)?.discarded); return { name: `donor__${record.id}`, sequence: record.sequence.replaceAll("-", "").toUpperCase().replaceAll("U", "T"), color: contaminant ? CONTAMINATION_TIP_LEGEND[1].color : CONTAMINATION_TIP_LEGEND[2].color, category: contaminant ? CONTAMINATION_TIP_LEGEND[1].label : CONTAMINATION_TIP_LEGEND[2].label }; }),
      ];
      const aligned = await runScalableMsa(rows.map((row) => row.sequence), runAlivibeMsa, undefined, 3, "nucleotide");
      const fasta = rows.map((row, index) => `>${row.name}\n${aligned[index]}\n`).join(""), newick = await runFastTree(fasta);
      setContaminationPhylogeny({ sample: targetSample, fasta, newick,
        leafMetadata: Object.fromEntries(rows.map((row) => [row.name, { familyCount: 1, color: row.color, category: row.category }])) });
    } catch (cause) { setContaminationPhylogenyError(cause instanceof Error ? cause.message : String(cause)); }
    finally { setContaminationPhylogenyBusy(false); }
  }

  async function returnFromAlivibe(session: AlivibeSession) {
    if (alivibeRef.current?.token !== session.token) throw new Error("This Alivibe editor belongs to an expired round trip.");
    if (session.popup.closed) throw new Error("The Alivibe editor has been closed.");
    const current = effectiveAlignment(bundleRef.current, session.sample, session.variant).fasta;
    if (!current) throw new Error("The originating nucleotide alignment is no longer loaded.");
    assertAlivibeRoundTripTarget({ groupKey: `${session.sample}/${session.variant}`, alignmentFingerprint: session.baselineFingerprint }, { groupKey: `${sampleRef.current}/${session.variant}`, alignmentFingerprint: inspectAlignment(current, 1).fingerprint });
    const returned = readAlivibeNucleotideFasta(session.popup);
    await installCorrection(session.variant, session.sample, session.baseline, returned.fasta, returned.frameOffset, `Alivibe ${returned.sourceRevision.slice(0, 12)}`);
    alivibeRef.current = null; window.focus(); session.popup.close();
  }

  function openAlivibe(variant: AlignmentVariant) {
    const existing = alivibeRef.current; if (existing && !existing.popup.closed) { existing.popup.focus(); setAlignmentError("An Alivibe round trip is already open."); return; }
    const targetSample = sampleRef.current, current = effectiveAlignment(bundleRef.current, targetSample, variant); if (!current.fasta) return;
    setAlignmentError(""); setAlignmentStatus("Opening the bundled Alivibe editor…");
    const applicationBase = new URL(import.meta.env.BASE_URL, document.baseURI), editorUrl = new URL("tools/alivibe.html", applicationBase);
    editorUrl.searchParams.set("swigBridge", String(ALIVIBE_BRIDGE_VERSION)); editorUrl.searchParams.set("source", ALIVIBE_SOURCE_REVISION.slice(0, 12)); editorUrl.searchParams.set("release", "0.3.10");
    const token = globalThis.crypto?.randomUUID?.() ?? `${Date.now()}-${Math.random().toString(36).slice(2)}`;
    const popup = window.open(editorUrl.href, `webporpid-alivibe-${token}`, "popup,width=1500,height=920") as AlivibeEditorWindow | null;
    if (!popup) { setAlignmentError("The browser blocked the Alivibe window. Allow pop-ups for this site and try again."); return; }
    const session: AlivibeSession = { token, popup, sample: targetSample, variant, baseline: current.fasta, baselineFingerprint: inspectAlignment(current.fasta, 1).fingerprint, frameOffset: current.frameOffset };
    alivibeRef.current = session; let attempts = 0;
    const initialize = () => {
      attempts += 1;
      if (popup.closed) { if (session.timer) window.clearInterval(session.timer); if (alivibeRef.current?.token === token) alivibeRef.current = null; setAlignmentStatus("Alivibe closed without returning an alignment."); return; }
      try {
        const bridge = getAlivibeBridge(popup); if (!bridge) { if (attempts < 240) return; throw new Error("The bundled Alivibe bridge did not become ready."); }
        bridge.installMsaRunner(createMsaJob); bridge.installFastTreeRunner((fasta, mode) => runFastTree(fasta, mode));
        const loaded = loadAlivibeNucleotideFasta(popup, session.baseline, session.frameOffset); assertAlivibeInitialLoad(session.baseline, loaded);
        const controls = popup.document.getElementById("controls"); if (!controls) throw new Error("Alivibe loaded without its controls.");
        if (!popup.document.getElementById("webporpid-return-control")) {
          const group = popup.document.createElement("div"); group.id = "webporpid-return-control"; group.className = "control-group";
          const label = popup.document.createElement("label"); label.textContent = "webPORPID round trip";
          const button = popup.document.createElement("button"); button.type = "button"; button.textContent = "Return alignment to webPORPID"; button.className = "active";
          button.onclick = () => void returnFromAlivibe(session).catch((cause) => { setAlignmentError(cause instanceof Error ? cause.message : String(cause)); window.focus(); });
          group.append(label, button); controls.prepend(group);
        }
        setAlignmentStatus(`Nucleotide alignment loaded in bundled Alivibe. Edit it, then press Return alignment to webPORPID. The original alignment will remain preserved.`);
        if (session.timer) window.clearInterval(session.timer);
      } catch (cause) {
        setAlignmentError(cause instanceof Error ? cause.message : String(cause)); setAlignmentStatus("Alivibe round trip stopped before editing.");
        if (alivibeRef.current?.token === token) alivibeRef.current = null; if (session.timer) window.clearInterval(session.timer);
      }
    };
    popup.addEventListener("load", initialize, { once: true }); session.timer = window.setInterval(initialize, 250);
  }

  async function importCorrected(variant: AlignmentVariant, file?: File) {
    const current = effectiveAlignment(bundleRef.current, sampleRef.current, variant); if (!file || !current.fasta) return;
    await installCorrection(variant, sampleRef.current, current.fasta, await file.text(), current.frameOffset, `Imported corrected FASTA · ${file.name}`);
  }

  function alignmentBlock(variant: AlignmentVariant) {
    const current = variant === "collapsed" ? collapsed : variant === "uncollapsed" ? uncollapsed : functional;
    const tree = variant === "collapsed" ? collapsedTree : variant === "uncollapsed" ? uncollapsedTree : functionalTree;
    const metadata = variant === "collapsed" ? collapsedTips : variant === "uncollapsed" ? uncollapsedTips : functionalTips;
    const reference = variant === "functional" ? functionalRefSequence : refSequence, busy = alignmentBusy === `${sample}/${variant}`;
    const kicker = variant === "collapsed" ? "Default phylogeny" : variant === "uncollapsed" ? "Optional family-level phylogeny" : "Functional-sequence phylogeny";
    const heading = variant === "collapsed" ? "Collapsed haplotypes" : variant === "uncollapsed" ? "Uncollapsed consensus sequences" : "Functional collapsed variants";
    const description = variant === "collapsed"
      ? "Identical retained consensuses are collapsed; bubble area represents UMI-family count, never read count. Minimum agreement remains a per-family property and is available in the uncollapsed view."
      : variant === "uncollapsed" ? "Every retained UMI-family consensus is shown separately and can be colored by its family minimum agreement."
        : "Every collapsed variant that passed the functional filter retains its abundance-ranked collapsed name and is shown in the joint codon-aware alignment, clipped exactly to the first and last nucleotide of the functional reference.";
    if (!current.fasta) return <div className="empty-state">No {variant} nucleotide alignment is available for this sample.</div>;
    const stale = Boolean(current.edit?.treeStale || (current.edit?.treeFingerprint && current.edit.treeFingerprint !== current.edit.editedFingerprint));
    return <section className="phylogeny-block" key={variant}><header><div><span className="section-kicker">{kicker}</span><h3>{heading}</h3><p>{description}</p><MethodLink topic={variant === "functional" ? "functional" : variant === "collapsed" ? "collapse" : "phylogeny"} /></div><div><button type="button" className={tree ? "" : "primary"} disabled={busy} onClick={() => void inferTree(variant)}>{busy ? "Running FastTree…" : tree ? "Recalculate tree" : "Infer tree"}</button></div></header>
      <div className="alignment-edit-bar"><button type="button" disabled={busy} onClick={() => openAlivibe(variant)}>Open in Alivibe ↗</button><label className="button-like">Import edited FASTA<input type="file" accept=".fasta,.fa,.fas,.fna,text/plain" onChange={(event) => { const file = event.target.files?.[0]; event.currentTarget.value = ""; void importCorrected(variant, file); }} /></label><label><span>Protein frame</span><select value={current.frameOffset} onChange={(event) => setFrameOffset(variant, Number(event.target.value) as AlignmentFrameOffset)}><option value="0">Start at nucleotide column 1</option><option value="1">Start at nucleotide column 2</option><option value="2">Start at nucleotide column 3</option></select></label>{current.edit && <button type="button" onClick={() => resetCorrection(variant)}>Discard edit</button>}<span>{current.edit ? `Edited copy · ${current.edit.source}` : "Pipeline alignment · original preserved"}</span></div>
      {current.edit && <div className="alignment-edited-notice"><strong>Alignment edited</strong><span>This separate copy is stored in the session file. The pipeline-generated alignment has not been overwritten.{stale ? " Recalculate the tree when the alignment is ready." : ""}</span></div>}
      <AlignmentTreeViewer fasta={current.fasta} newick={tree} alphabet={alphabet} onAlphabetChange={setAlphabet} frameOffset={current.frameOffset} referenceSequence={reference} leafMetadata={metadata} collapsed={variant !== "uncollapsed"} treeStale={stale} name={`${safeDatasetName(bundle.config.dataset)}-${safeDatasetName(sample)}-${variant}`} />
    </section>;
  }

  return <main className="results-page">
    <section className="results-hero"><div><span className="section-kicker">Loaded analysis</span><h1>{bundle.config.dataset}</h1><p>{bundle.provenance.inputName} · {bundle.provenance.createdUtc} · {bundle.provenance.workers} workers</p>{(!allSamples || donorView) && <button type="button" className="overview-return" onClick={() => { setDonorView(""); setAllSamples(true); }}>← Across-sample overview</button>}</div><div className="result-actions"><select aria-label="Sample, donor, or across-sample overview" value={donorView ? `__donor__${donorView}` : allSamples ? "__all_samples__" : sample} onChange={(event) => { const value = event.target.value; if (value === "__all_samples__") { setDonorView(""); setAllSamples(true); } else if (value.startsWith("__donor__")) { setDonorView(value.slice(9)); setAllSamples(false); } else { setDonorView(""); setSample(value); setAllSamples(false); setTab("overview"); } }}><option value="__all_samples__">Across-sample overview</option>{donorIds.length > 0 && <optgroup label="Donors">{donorIds.map((donor) => <option value={`__donor__${donor}`} key={donor}>Donor · {donor}</option>)}</optgroup>}<optgroup label="Samples">{bundle.summaries.map((row) => <option value={row.sample} key={row.sample}>{row.sample}</option>)}</optgroup></select><ExportMenu bundle={bundle} sample={sample} allOnly={allSamples || Boolean(donorView)} /><button className="primary" type="button" onClick={onSaveResults}>Save results file</button></div></section>
    <OptionalStagePanel bundle={bundle} progress={deferredProgress} onCompute={(stage) => void computeOptionalStage(stage)} onSkip={skipDeferredStage} />
    {deferredError && <div className="error-box" role="alert">Optional-stage computation failed: {deferredError}</div>}
    {allSamples && !donorView && <AllSampleOverview bundle={bundle} onOpenSample={(selected) => { setDonorView(""); setSample(selected); setAllSamples(false); setTab("overview"); }} />}
    {donorView && <DonorView bundle={bundle} donorId={donorView} />}
    {sampleMode && <nav className="result-tabs">{(["overview", "families", "sequences", "contamination", "alignment", "log"] as Tab[]).map((value) => <button type="button" className={tab === value ? "active" : ""} onClick={() => setTab(value)} key={value}>{value}</button>)}</nav>}
    {sampleMode && <>{tab === "overview" && <section className="result-section">
      <div className="metric-grid"><article><span>Demultiplexed reads</span><strong>{summary?.demultiplexedReads.toLocaleString() ?? "0"}</strong></article><article><span>Selected reads</span><strong>{currentOverview?.selectedReads.toLocaleString() ?? "0"}</strong></article><article><span>Subsampled reads</span><strong>{currentOverview ? `${currentOverview.downsampledReads.toLocaleString()} · ${pct(currentOverview.downsampledPercent)}` : "—"}</strong></article><article><span>Observed UMI families</span><strong>{summary?.observedUmis.toLocaleString() ?? "0"}</strong></article><article><span>Consensus families</span><strong>{summary?.consensusSequences.toLocaleString() ?? "0"}</strong></article><article><span>Collapsed haplotypes</span><strong>{summary?.collapsedSequences?.toLocaleString() ?? bundle.collapseGroups?.[sample]?.length.toLocaleString() ?? "—"}</strong></article><article><span>Functional collapsed variants</span><strong>{summary?.functionalPassed?.toLocaleString() ?? "—"}</strong></article><article><span>Postproc passed families</span><strong>{summary?.postprocPassed?.toLocaleString() ?? "—"}</strong></article><article><span>Artefact cutoff</span><strong>{summary?.artefactCutoff ?? "—"}</strong></article></div>
      {postprocessingDone && <Filters bundle={bundle} sample={sample} />}
      <div className="statistics-grid"><DualStatsTable title="PORPID call statistics" description="Each observed UMI is counted once. The aggregate BPB-reject bucket contributes reads / CCS but not an artificial UMI count; family sizes expose heteroduplex and LDA impact at both molecule and read level." rows={porpidStats} />{postprocessingDone && <DualStatsTable title="Consensus and post-processing filters" description="Consensus-filter rows report every matching rejection category; a family can match more than one row." rows={postprocStats} familyLabel="Consensus families" />}{collapseDone && <DualStatsTable title="Functional-filter results" description="Functional filtering is evaluated once per collapsed variant. Variant counts are unweighted; represented-family counts preserve the number of UMI families in each variant and never use read counts." rows={functionalStats} familyLabel="Collapsed variants" readLabel="Represented UMI families" />}</div>
      {postprocessingDone ? <div className="chart-grid report-figures"><article><header><h3>UMI family size × UMI length</h3><p>Family-size, UMI-length, probabilistic class, quantile, and artefact decisions.</p></header><UmiDecisionPlot families={families} artefactCutoff={summary?.artefactCutoff ?? 0} agreementThreshold={sampleConfig?.agreementOverride ?? bundle.config.parameters.agreementThreshold} outlierQuantile={sampleConfig?.outlierQuantileOverride ?? bundle.config.parameters.outlierQuantile} /></article><article><header><h3>Artefact-cutoff decision</h3><p>Deterministic jitter below the non-outlier quantile with threshold guides.</p></header><ArtefactDecisionPlot families={families} artefactCutoff={summary?.artefactCutoff ?? 0} agreementThreshold={sampleConfig?.agreementOverride ?? bundle.config.parameters.agreementThreshold} outlierQuantile={sampleConfig?.outlierQuantileOverride ?? bundle.config.parameters.outlierQuantile} artefactFraction={sampleConfig?.artefactFractionOverride ?? bundle.config.parameters.artefactFraction} /></article><article><header><h3>Low-agreement positions</h3><p>One longest-run minimum site per non-artefact family, sized by homopolymer run and colored by modal base.</p></header><AgreementPositionPlot consensuses={consensuses} threshold={sampleConfig?.agreementOverride ?? bundle.config.parameters.agreementThreshold} minimumFamilySize={summary?.artefactCutoff ?? 0} /></article><article><header><h3>P(APOBEC|mutations) · classical MDS</h3><p>Pairwise non-gap Hamming distances; point area follows family size and color follows the stored APOBEC posterior.</p></header><MdsApobecPlot records={sampleRecords} /></article><article><header><h3>UMI dinucleotide frequencies</h3><p>Family-unweighted frequencies beside frequencies weighted by reads within each UMI family.</p></header><DinucleotideHeatmaps families={families} /></article></div> : <div className="empty-state">Downstream decision plots have not been computed. Use a Compute through action above to generate them from the stored consensus calls.</div>}
      <TimingSummary bundle={bundle} /><article className="provenance-card"><h3>Audit trail</h3><dl><div><dt>Input SHA-256</dt><dd>{bundle.provenance.inputSha256}</dd></div><div><dt>Config SHA-256</dt><dd>{bundle.provenance.configSha256}</dd></div><div><dt>Engine</dt><dd>{bundle.provenance.engine}</dd></div><div><dt>Processing source</dt><dd>{bundle.provenance.upstreamBranch}@{bundle.provenance.upstreamCommit.slice(0, 12)}</dd></div>{bundle.inputMappings?.map((mapping) => <div key={`${mapping.slot}-${mapping.uploadedName}`}><dt>{mapping.role} · {mapping.slot}</dt><dd>{mapping.uploadedName}{mapping.expectedName && mapping.expectedName !== mapping.uploadedName ? ` → ${mapping.expectedName}` : ""}</dd></div>)}</dl><h4>Interactive threshold decisions</h4><ThresholdAuditTrail bundle={bundle} sample={sample} /><h4>Alignment and phylogeny changes</h4><AlignmentAuditTable bundle={bundle} sample={sample} /></article>
    </section>}
    {tab === "families" && <section className="result-section"><div className="table-heading"><div><h2>UMI family decisions</h2><p>Probabilistic offspring assignment, heteroduplex check, length and family-size gates.</p><MethodLink topic="umi" /></div><input placeholder="Search UMI, parent, or class" value={familyQuery} onChange={(event) => setFamilyQuery(event.target.value)} /></div><div className="table-scroll"><table><thead><tr><SortHeader label="UMI" column="umi" state={familySort} onChange={setFamilySort} /><SortHeader label="fs" column="familySize" state={familySort} onChange={setFamilySort} /><SortHeader label="classification" column="classification" state={familySort} onChange={setFamilySort} /><SortHeader label="parent" column="mostLikelyParent" state={familySort} onChange={setFamilySort} /><SortHeader label="posterior" column="posteriorProbability" state={familySort} onChange={setFamilySort} /><SortHeader label="minag" column="minimumAgreement" state={familySort} onChange={setFamilySort} /></tr></thead><tbody>{page(filteredFamilies, familyPage).map((row) => <tr key={`${row.sampleIndex}-${row.umi}`}><td><code>{row.umi}</code></td><td>{row.familySize}</td><td><span className={`status ${row.disposition === "likely_real" ? "pass" : "reject"}`}>{row.disposition}</span></td><td><code>{row.mostLikelyParent}</code></td><td>{row.posteriorProbability.toFixed(6)}</td><td>{row.minimumAgreement?.toFixed(2) ?? "—"}</td></tr>)}</tbody></table></div><Pager page={familyPage} count={filteredFamilies.length} onChange={setFamilyPage} /></section>}
    {tab === "sequences" && <section className="result-section"><div className="table-heading"><div><h2>Consensus and post-processing</h2><p>Every filter decision remains inspectable; FASTA exports use the exact stored sequences.</p><MethodLink topic="consensus" /></div><input placeholder="Search ID or UMI" value={query} onChange={(event) => setQuery(event.target.value)} /></div>{!postprocessingDone && <div className="viewer-status warning">Downstream filtering has not been computed; no filter decisions are being presented as passes or rejects.</div>}<div className="table-scroll"><table><thead><tr><SortHeader label="Sequence" column="id" state={sequenceSort} onChange={setSequenceSort} /><SortHeader label="fs" column="familySize" state={sequenceSort} onChange={setSequenceSort} /><SortHeader label="minag" column="minimumAgreement" state={sequenceSort} onChange={setSequenceSort} /><SortHeader label="panel score" column="panelScore" state={sequenceSort} onChange={setSequenceSort} /><SortHeader label="filters" column="filters" state={sequenceSort} onChange={setSequenceSort} /><SortHeader label="APOBEC p" column="apobec" state={sequenceSort} onChange={setSequenceSort} /><SortHeader label="reason" column="reason" state={sequenceSort} onChange={setSequenceSort} /></tr></thead><tbody>{page(records, sequencePage).map((row) => <tr key={row.id}><td><code title={row.id}>{row.id}</code></td><td>{row.familySize}</td><td>{row.minimumAgreement.toFixed(2)}</td><td>{row.panelScore.toFixed(2)}</td><td><span className={`status ${row.artefactPass && row.agreementPass && row.contaminationPass && row.panelPass ? "pass" : "reject"}`}>{row.artefactPass && row.agreementPass && row.contaminationPass && row.panelPass ? "pass" : "reject"}</span></td><td>{row.apobec?.posteriorGaInflated.toFixed(3) ?? "—"}</td><td>{row.rejectionReasons.join("; ") || "—"}</td></tr>)}</tbody></table></div><Pager page={sequencePage} count={records.length} onChange={setSequencePage} /></section>}
    {tab === "contamination" && <section className="result-section"><div className="table-heading"><div><h2>Contamination report</h2><p>One decision is shown per consensus family. A primary decision takes precedence; the wider suspect pass appears only when no primary call exists.</p><MethodLink topic="contamination" /></div><span>{contaminationRows.length.toLocaleString()} unique families</span></div>{!contaminationDone && <div className="viewer-status warning">Contamination checks have not been computed. An empty table here is not a pass decision. Other downstream outputs remain valid but explicitly unfiltered at this gate.</div>}<div className="table-scroll"><table><thead><tr><SortHeader label="Sequence" column="sequenceId" state={contaminationSort} onChange={setContaminationSort} /><SortHeader label="Nearest non-self" column="nearestNonselfVariant" state={contaminationSort} onChange={setContaminationSort} /><SortHeader label="distance" column="nearestNonselfDistance" state={contaminationSort} onChange={setContaminationSort} /><SortHeader label="decision" column="decision" state={contaminationSort} onChange={setContaminationSort} /></tr></thead><tbody>{page(contaminationRows, contaminationPage).map((row) => <tr key={row.sequenceId}><td><code>{row.sequenceId}</code></td><td>{row.nearestNonselfVariant}</td><td>{row.nearestNonselfDistance.toPrecision(5)}</td><td><span className={`status ${row.discarded ? "reject" : "warn"}`}>{row.discarded ? "discarded contaminant" : row.suspectOnly ? "suspect pass" : "reported, retained"}</span></td></tr>)}</tbody></table></div><Pager page={contaminationPage} count={contaminationRows.length} onChange={setContaminationPage} />
      <section className="contamination-phylogeny"><header><div><span className="section-kicker">On-demand cross-check</span><h3>Contamination references + donor sequences</h3><p>Build a joint nucleotide alignment and double-precision FastTree phylogeny. Tips are colored as panel references, donor contaminants discarded by the filter, or retained donor sequences.</p></div><button type="button" disabled={contaminationPhylogenyBusy || !contaminationDone} onClick={() => void buildContaminationPhylogeny()}>{!contaminationDone ? "Compute contamination checks first" : contaminationPhylogenyBusy ? "Aligning and inferring…" : contaminationPhylogeny ? "Rebuild tree + alignment" : "Build tree + alignment"}</button></header>
        {contaminationPhylogenyError && <div className="error-box" role="alert">{contaminationPhylogenyError}</div>}
        {contaminationPhylogeny?.sample === sample && <AlignmentTreeViewer fasta={contaminationPhylogeny.fasta} newick={contaminationPhylogeny.newick} alphabet={alphabet} onAlphabetChange={setAlphabet} referenceSequence={parseFasta(contaminationPhylogeny.fasta)[0]?.sequence} leafMetadata={contaminationPhylogeny.leafMetadata} tipLegend={CONTAMINATION_TIP_LEGEND} variantLabel="Contamination-reference phylogeny" name={`${safeDatasetName(bundle.config.dataset)}-${safeDatasetName(sample)}-contamination`} />}
      </section>
    </section>}
    {tab === "alignment" && <section className="result-section"><div className="alignment-switch"><div><h2>Tree + alignment workbench</h2><p>Collapsed haplotypes are the default phylogeny. Amino acid is translated directly from each active nucleotide alignment.</p><MethodLink topic="phylogeny" /></div><div><button className={alphabet === "nt" ? "active" : ""} onClick={() => setAlphabet("nt")}>Nucleotide alignment</button><button className={alphabet === "aa" ? "active" : ""} onClick={() => setAlphabet("aa")}>Amino-acid alignment</button></div></div>
      {edited && <div className="alignment-edited-notice"><strong>Edited alignments are present</strong><span>They are stored as separate session copies; original pipeline alignments remain available.</span></div>}{alignmentStatus && <div className="viewer-status">{alignmentStatus}</div>}{alignmentError && <div className="error-box" role="alert">{alignmentError}</div>}
      {alignmentBlock("collapsed")}
      {uncollapsed.fasta && (!showUncollapsed ? <button type="button" className="show-uncollapsed" onClick={() => setShowUncollapsed(true)}>Open uncollapsed family-level phylogeny + alignment</button> : <>{alignmentBlock("uncollapsed")}<button type="button" onClick={() => setShowUncollapsed(false)}>Close uncollapsed view</button></>)}
      {functional.fasta && (!showFunctional ? <button type="button" className="show-uncollapsed" onClick={() => setShowFunctional(true)}>Open all functional sequences + infer phylogeny</button> : <>{alignmentBlock("functional")}<button type="button" onClick={() => setShowFunctional(false)}>Close functional-sequence view</button></>)}
      {collapseDone && sampleConfig?.functionalReference && !functional.fasta && <div className="empty-state">No collapsed variants passed the functional filter for this sample, so there is no functional-sequence alignment to display.</div>}
    </section>}
    {tab === "log" && <section className="result-section"><div className="table-heading"><div><h2>Run log</h2><p>Persistent stage summaries, input-slot mappings, interactive tree runs, and fallbacks stored inside the results file.</p></div></div><pre className="run-log">{bundle.log.join("\n")}</pre></section>}</>}
  </main>;
}
