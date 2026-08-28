import { useEffect, useMemo, useRef, useState } from "react";
import {
  ALIVIBE_BRIDGE_VERSION, ALIVIBE_SOURCE_REVISION, assertAlivibeInitialLoad, assertAlivibeRoundTripTarget,
  getAlivibeBridge, loadAlivibeNucleotideFasta, readAlivibeNucleotideFasta, type AlivibeEditorWindow, type AlivibeMsaJob,
} from "../alivibe-roundtrip";
import { runAlivibeMsa } from "../alivibe-msa-runtime";
import { alignmentKey, effectiveAlignment, inspectAlignment, validateCorrectedAlignment, type AlignmentFrameOffset, type AlignmentVariant } from "../alignment-utils";
import { runFastTree } from "../biowasm";
import { parseFasta } from "../config";
import { buildExportArchive } from "../export-archive";
import { exportComponent, type ExportKind, safeDatasetName } from "../result-file";
import type { CollapseGroup, ContaminationCall, PostprocRecord, ResultBundle, UmiFamily } from "../types";
import { AgreementPositionPlot, ArtefactDecisionPlot, DinucleotideHeatmaps, MdsApobecPlot, UmiDecisionPlot } from "./charts";
import { AlignmentTreeViewer, type LeafMetadata } from "./alignment-tree-viewer";

type Tab = "overview" | "families" | "sequences" | "contamination" | "alignment" | "log";
type SortDirection = "asc" | "desc";
interface SortState<K extends string> { key: K; direction: SortDirection }
const PAGE_SIZE = 250;

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

function ExportMenu({ bundle, sample }: { bundle: ResultBundle; sample: string }) {
  const [kind, setKind] = useState<ExportKind>("consensus-fasta"), [exportingAll, setExportingAll] = useState(false);
  const labels: Array<[ExportKind, string]> = [["consensus-fasta", "Consensus FASTA"], ["passed-consensus-fasta", "Passed consensus FASTA"],
    ["rejected-consensus-fasta", "Rejected consensus FASTA"], ["trimmed-nt-fasta", "Trimmed nucleotide FASTA"], ["trimmed-aa-fasta", "Trimmed amino-acid FASTA"],
    ["family-csv", "UMI family CSV"], ["low-agreement-csv", "Low-agreement CSV"], ["contamination-csv", "Contamination CSV"],
    ["postproc-csv", "Postproc CSV"], ["apobec-csv", "APOBEC CSV"], ["collapse-csv", "Collapse membership CSV"],
    ["nucleotide-alignment", "Collapsed nucleotide alignment"], ["protein-alignment", "Collapsed protein alignment"], ["newick", "Collapsed Newick tree"],
    ["uncollapsed-nucleotide-alignment", "Uncollapsed nucleotide alignment"], ["uncollapsed-protein-alignment", "Uncollapsed protein alignment"], ["uncollapsed-newick", "Uncollapsed Newick tree"], ["log", "Run log"]];
  return <div className="export-menu"><select value={kind} onChange={(event) => setKind(event.target.value as ExportKind)}>{labels.map(([value, label]) => <option value={value} key={value}>{label}</option>)}</select>
    <button type="button" onClick={() => { const result = exportComponent(bundle, kind, sample); downloadData(`${safeDatasetName(bundle.config.dataset)}-${safeDatasetName(sample)}.${result.extension}`, result.text, result.mime); }}>Export</button>
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

function editWarnings(inspection: ReturnType<typeof validateCorrectedAlignment>): string[] {
  const warnings: string[] = [];
  if (inspection.removedRows.length) warnings.push(`${inspection.removedRows.length} sequence rows deleted`);
  if (inspection.addedRows.length) warnings.push(`${inspection.addedRows.length} sequence rows added or renamed`);
  if (inspection.substitutedNucleotides) warnings.push(`${inspection.substitutedNucleotides} base positions substituted`);
  if (inspection.insertedNucleotides) warnings.push(`${inspection.insertedNucleotides} nucleotide characters inserted`);
  if (inspection.removedNucleotides) warnings.push(`${inspection.removedNucleotides} nucleotide characters deleted`);
  return warnings;
}

function referenceSequence(bundle: ResultBundle, sample: string): string | undefined {
  const fasta = bundle.referenceAlignments?.[alignmentKey(sample)];
  if (!fasta) return undefined;
  try { return parseFasta(fasta)[0]?.sequence; } catch { return undefined; }
}

function collapsedMetadata(groups: CollapseGroup[] | undefined): Record<string, LeafMetadata> {
  return Object.fromEntries((groups ?? []).map((group) => [group.representativeId, { familyCount: group.familyCount }]));
}

function uncollapsedMetadata(records: PostprocRecord[]): Record<string, LeafMetadata> {
  return Object.fromEntries(records.filter((record) => record.alignedNt).map((record) => [record.id, { familyCount: 1, minimumAgreement: record.minimumAgreement }]));
}

export function ResultsExplorer({ bundle, onSaveResults, onBundleChange }: { bundle: ResultBundle; onSaveResults(): void; onBundleChange(bundle: ResultBundle): void }) {
  const [tab, setTab] = useState<Tab>("overview"), [sample, setSample] = useState(bundle.summaries[0]?.sample ?? bundle.config.samples[0]?.name ?? "");
  const [query, setQuery] = useState(""), [familyQuery, setFamilyQuery] = useState(""), [alphabet, setAlphabet] = useState<"nt" | "aa">("nt"), [showUncollapsed, setShowUncollapsed] = useState(false);
  const [familyPage, setFamilyPage] = useState(0), [sequencePage, setSequencePage] = useState(0), [contaminationPage, setContaminationPage] = useState(0);
  const [familySort, setFamilySort] = useState<SortState<string>>({ key: "familySize", direction: "desc" });
  const [sequenceSort, setSequenceSort] = useState<SortState<string>>({ key: "minimumAgreement", direction: "asc" });
  const [contaminationSort, setContaminationSort] = useState<SortState<string>>({ key: "nearestNonselfDistance", direction: "asc" });
  const [alignmentStatus, setAlignmentStatus] = useState(""), [alignmentError, setAlignmentError] = useState(""), [alignmentBusy, setAlignmentBusy] = useState("");
  const bundleRef = useRef(bundle), sampleRef = useRef(sample), changeRef = useRef(onBundleChange), alivibeRef = useRef<AlivibeSession | null>(null);
  bundleRef.current = bundle; sampleRef.current = sample; changeRef.current = onBundleChange;

  const summary = bundle.summaries.find((row) => row.sample === sample), sampleConfig = bundle.config.samples.find((row) => row.name === sample);
  const families = useMemo(() => bundle.umiFamilies.filter((row) => row.sample === sample), [bundle, sample]);
  const consensuses = useMemo(() => bundle.consensuses.filter((row) => row.sample === sample), [bundle, sample]);
  const filteredFamilies = useMemo(() => sorted(families.filter((row) => !familyQuery || row.umi.includes(familyQuery.toUpperCase()) || row.mostLikelyParent.includes(familyQuery.toUpperCase()) || row.disposition.toLowerCase().includes(familyQuery.toLowerCase())), familySort,
    (row: UmiFamily, key) => key === "classification" ? row.disposition : row[key as keyof UmiFamily]), [families, familyQuery, familySort]);
  const sampleRecords = useMemo(() => bundle.records.filter((row) => row.sample === sample), [bundle, sample]);
  const records = useMemo(() => sorted(sampleRecords.filter((row) => !query || row.id.toLowerCase().includes(query.toLowerCase()) || row.umi.includes(query.toUpperCase())), sequenceSort,
    (row: PostprocRecord, key) => key === "filters" ? Number(row.artefactPass && row.agreementPass && row.contaminationPass && row.panelPass) : key === "apobec" ? row.apobec?.posteriorGaInflated : key === "reason" ? row.rejectionReasons.join(";") : row[key as keyof PostprocRecord]), [sampleRecords, query, sequenceSort]);
  const contaminationRows = useMemo(() => sorted(bundle.contamination.filter((row) => row.sample === sample), contaminationSort,
    (row: ContaminationCall, key) => key === "decision" ? row.discarded ? 2 : row.suspectOnly ? 1 : 0 : row[key as keyof ContaminationCall]), [bundle, sample, contaminationSort]);
  const collapsed = useMemo(() => effectiveAlignment(bundle, sample, "collapsed"), [bundle, sample]);
  const uncollapsed = useMemo(() => effectiveAlignment(bundle, sample, "uncollapsed"), [bundle, sample]);
  const collapsedTree = collapsed.edit?.treeNewick ?? bundle.trees[collapsed.key], uncollapsedTree = uncollapsed.edit?.treeNewick ?? bundle.trees[uncollapsed.key];
  const refSequence = useMemo(() => referenceSequence(bundle, sample), [bundle, sample]);
  const collapsedTips = useMemo(() => collapsedMetadata(bundle.collapseGroups?.[sample]), [bundle, sample]);
  const uncollapsedTips = useMemo(() => uncollapsedMetadata(sampleRecords), [sampleRecords]);
  const edited = Boolean(collapsed.edit || uncollapsed.edit);

  useEffect(() => { setFamilyPage(0); setSequencePage(0); setContaminationPage(0); setAlignmentStatus(""); setAlignmentError(""); setShowUncollapsed(false); }, [sample]);
  useEffect(() => setFamilyPage(0), [familyQuery, familySort]); useEffect(() => setSequencePage(0), [query, sequenceSort]); useEffect(() => setContaminationPage(0), [contaminationSort]);
  useEffect(() => () => { const session = alivibeRef.current; if (session?.timer) window.clearInterval(session.timer); if (session && !session.popup.closed) session.popup.close(); }, []);
  const page = <T,>(rows: T[], index: number) => rows.slice(index * PAGE_SIZE, (index + 1) * PAGE_SIZE);

  async function installCorrection(variant: AlignmentVariant, targetSample: string, baseline: string, corrected: string, frameOffset: AlignmentFrameOffset, source: string) {
    setAlignmentBusy(`${targetSample}/${variant}`); setAlignmentError("");
    try {
      const inspected = validateCorrectedAlignment(baseline, corrected), warnings = editWarnings(inspected);
      if (warnings.length && !window.confirm(`This alignment edit changes biological content:\n\n• ${warnings.join("\n• ")}\n\nThe pipeline-generated alignment will remain preserved. Save this separate edit?`)) {
        setAlignmentStatus("Alignment edit cancelled; the stored alignment was not changed."); return;
      }
      const currentBundle = bundleRef.current, key = alignmentKey(targetSample, variant), original = currentBundle.alignments[key] ?? (variant === "uncollapsed" ? currentBundle.alignments[alignmentKey(targetSample)] : undefined);
      if (!original) throw new Error("The original nucleotide alignment is unavailable.");
      const current = effectiveAlignment(currentBundle, targetSample, variant), existingTree = current.edit?.treeNewick ?? currentBundle.trees[key];
      const next: ResultBundle = { ...currentBundle, alignmentEdits: { ...currentBundle.alignmentEdits, [key]: {
        fasta: inspected.fasta, frameOffset, baselineFingerprint: inspectAlignment(original, 1).fingerprint, editedFingerprint: inspected.fingerprint,
        source, savedUtc: new Date().toISOString(), treeNewick: existingTree,
        treeFingerprint: current.edit?.treeFingerprint ?? (existingTree ? inspectAlignment(current.fasta ?? original, 1).fingerprint : undefined),
        treeStale: Boolean(existingTree) && inspected.fingerprint !== (current.edit?.treeFingerprint ?? inspectAlignment(current.fasta ?? original, 1).fingerprint), warnings,
      } } };
      bundleRef.current = next; changeRef.current(next);
      setAlignmentStatus(`Saved a separate ${variant} alignment edit: ${inspected.rows.toLocaleString()} rows × ${inspected.columns.toLocaleString()} columns.${warnings.length ? ` Warning: ${warnings.join("; ")}.` : " Gap placement changed; nucleotide content is unchanged."} Recalculate the tree when ready.`);
    } catch (cause) { setAlignmentError(cause instanceof Error ? cause.message : String(cause)); }
    finally { setAlignmentBusy(""); }
  }

  function setFrameOffset(variant: AlignmentVariant, frameOffset: AlignmentFrameOffset) {
    const currentBundle = bundleRef.current, targetSample = sampleRef.current, current = effectiveAlignment(currentBundle, targetSample, variant), key = current.key;
    if (!current.fasta) return;
    const fingerprint = inspectAlignment(current.fasta, 1).fingerprint, original = currentBundle.alignments[key] ?? current.fasta;
    const next: ResultBundle = { ...currentBundle, alignmentEdits: { ...currentBundle.alignmentEdits, [key]: {
      fasta: current.fasta, frameOffset, baselineFingerprint: inspectAlignment(original, 1).fingerprint, editedFingerprint: fingerprint,
      source: current.edit?.source ?? "Translation frame selection", savedUtc: new Date().toISOString(), treeNewick: current.edit?.treeNewick ?? currentBundle.trees[key],
      treeFingerprint: current.edit?.treeFingerprint ?? (currentBundle.trees[key] ? fingerprint : undefined), treeStale: current.edit?.treeStale ?? false, warnings: current.edit?.warnings,
    } } };
    bundleRef.current = next; changeRef.current(next); setAlignmentStatus(`Protein translation now starts at nucleotide column ${frameOffset + 1}; this choice is stored without changing the original alignment.`);
  }

  function resetCorrection(variant: AlignmentVariant) {
    const currentBundle = bundleRef.current, key = alignmentKey(sampleRef.current, variant), edits = { ...currentBundle.alignmentEdits };
    delete edits[key]; const next: ResultBundle = { ...currentBundle, alignmentEdits: Object.keys(edits).length ? edits : undefined };
    bundleRef.current = next; changeRef.current(next); setAlignmentError(""); setAlignmentStatus(`Restored the pipeline-generated ${variant} alignment and its stored tree.`);
  }

  async function inferTree(variant: AlignmentVariant) {
    const currentBundle = bundleRef.current, targetSample = sampleRef.current, current = effectiveAlignment(currentBundle, targetSample, variant), key = current.key;
    if (!current.fasta) return;
    setAlignmentBusy(`${targetSample}/${variant}`); setAlignmentError(""); setAlignmentStatus(`Inferring the ${variant} tree with bundled double-precision FastTree…`);
    try {
      const treeNewick = await runFastTree(current.fasta), fingerprint = inspectAlignment(current.fasta, 1).fingerprint;
      let next: ResultBundle;
      if (current.edit) next = { ...currentBundle, alignmentEdits: { ...currentBundle.alignmentEdits, [key]: { ...current.edit, treeNewick, treeFingerprint: fingerprint, treeStale: false, savedUtc: new Date().toISOString() } }, log: [...currentBundle.log, `${new Date().toISOString()} interactive FastTree: ${key}; edited alignment fingerprint ${fingerprint}`] };
      else next = { ...currentBundle, trees: { ...currentBundle.trees, [key]: treeNewick }, log: [...currentBundle.log, `${new Date().toISOString()} interactive FastTree: ${key}; alignment fingerprint ${fingerprint}`] };
      bundleRef.current = next; changeRef.current(next); setAlignmentStatus(`${variant === "collapsed" ? "Collapsed" : "Uncollapsed"} tree inference complete. Mutation mapping has been refreshed for this alignment.`);
    } catch (cause) { setAlignmentError(cause instanceof Error ? cause.message : String(cause)); }
    finally { setAlignmentBusy(""); }
  }

  function createMsaJob(sequences: string[], scoringMode: "nucleotide" | "amino-acid"): AlivibeMsaJob {
    const controller = new AbortController(); return { result: runAlivibeMsa(sequences, controller.signal, 3, scoringMode), cancel: () => controller.abort() };
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
    editorUrl.searchParams.set("swigBridge", String(ALIVIBE_BRIDGE_VERSION)); editorUrl.searchParams.set("source", ALIVIBE_SOURCE_REVISION.slice(0, 12)); editorUrl.searchParams.set("release", "0.3.4");
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
    const current = variant === "collapsed" ? collapsed : uncollapsed, tree = variant === "collapsed" ? collapsedTree : uncollapsedTree;
    const metadata = variant === "collapsed" ? collapsedTips : uncollapsedTips, busy = alignmentBusy === `${sample}/${variant}`;
    if (!current.fasta) return <div className="empty-state">No {variant} nucleotide alignment is available for this sample.</div>;
    const stale = Boolean(current.edit?.treeStale || (current.edit?.treeFingerprint && current.edit.treeFingerprint !== current.edit.editedFingerprint));
    return <section className="phylogeny-block" key={variant}><header><div><span className="section-kicker">{variant === "collapsed" ? "Default phylogeny" : "Optional family-level phylogeny"}</span><h3>{variant === "collapsed" ? "Collapsed haplotypes" : "Uncollapsed consensus sequences"}</h3><p>{variant === "collapsed" ? "Identical retained consensuses are collapsed; bubble area represents UMI-family count, never read count. Minimum agreement remains a per-family property and is available in the uncollapsed view." : "Every retained UMI-family consensus is shown separately and can be colored by its family minimum agreement."}</p></div><div><button type="button" className={tree ? "" : "primary"} disabled={busy} onClick={() => void inferTree(variant)}>{busy ? "Running FastTree…" : tree ? "Recalculate tree" : "Infer tree"}</button></div></header>
      <div className="alignment-edit-bar"><button type="button" disabled={busy} onClick={() => openAlivibe(variant)}>Open in Alivibe ↗</button><label className="button-like">Import edited FASTA<input type="file" accept=".fasta,.fa,.fas,.fna,text/plain" onChange={(event) => { const file = event.target.files?.[0]; event.currentTarget.value = ""; void importCorrected(variant, file); }} /></label><label><span>Protein frame</span><select value={current.frameOffset} onChange={(event) => setFrameOffset(variant, Number(event.target.value) as AlignmentFrameOffset)}><option value="0">Start at nucleotide column 1</option><option value="1">Start at nucleotide column 2</option><option value="2">Start at nucleotide column 3</option></select></label>{current.edit && <button type="button" onClick={() => resetCorrection(variant)}>Discard edit</button>}<span>{current.edit ? `Edited copy · ${current.edit.source}` : "Pipeline alignment · original preserved"}</span></div>
      {current.edit && <div className="alignment-edited-notice"><strong>Alignment edited</strong><span>This separate copy is stored in the session file. The pipeline-generated alignment has not been overwritten.{stale ? " Recalculate the tree when the alignment is ready." : ""}</span></div>}
      <AlignmentTreeViewer fasta={current.fasta} newick={tree} alphabet={alphabet} frameOffset={current.frameOffset} referenceSequence={refSequence} leafMetadata={metadata} collapsed={variant === "collapsed"} treeStale={stale} name={`${safeDatasetName(bundle.config.dataset)}-${safeDatasetName(sample)}-${variant}`} />
    </section>;
  }

  return <main className="results-page">
    <section className="results-hero"><div><span className="section-kicker">Loaded analysis</span><h1>{bundle.config.dataset}</h1><p>{bundle.provenance.inputName} · {bundle.provenance.createdUtc} · {bundle.provenance.workers} workers</p></div><div className="result-actions"><select aria-label="Sample" value={sample} onChange={(event) => setSample(event.target.value)}>{bundle.summaries.map((row) => <option key={row.sample}>{row.sample}</option>)}</select><ExportMenu bundle={bundle} sample={sample} /><button className="primary" type="button" onClick={onSaveResults}>Save results file</button></div></section>
    <nav className="result-tabs">{(["overview", "families", "sequences", "contamination", "alignment", "log"] as Tab[]).map((value) => <button type="button" className={tab === value ? "active" : ""} onClick={() => setTab(value)} key={value}>{value}</button>)}</nav>
    {tab === "overview" && <section className="result-section">
      <div className="metric-grid"><article><span>Demultiplexed</span><strong>{summary?.demultiplexedReads.toLocaleString() ?? "0"}</strong></article><article><span>Observed UMIs</span><strong>{summary?.observedUmis.toLocaleString() ?? "0"}</strong></article><article><span>Consensus</span><strong>{summary?.consensusSequences.toLocaleString() ?? "0"}</strong></article><article><span>Collapsed haplotypes</span><strong>{summary?.collapsedSequences?.toLocaleString() ?? bundle.collapseGroups?.[sample]?.length.toLocaleString() ?? "—"}</strong></article><article><span>Postproc passed</span><strong>{summary?.postprocPassed.toLocaleString() ?? "0"}</strong></article><article><span>Artefact cutoff</span><strong>{summary?.artefactCutoff ?? 0}</strong></article></div>
      <Filters bundle={bundle} sample={sample} />
      <div className="chart-grid report-figures"><article><header><h3>UMI family size × UMI length</h3><p>Family-size, UMI-length, probabilistic class, quantile, and artefact decisions.</p></header><UmiDecisionPlot families={families} artefactCutoff={summary?.artefactCutoff ?? 0} agreementThreshold={sampleConfig?.agreementOverride ?? bundle.config.parameters.agreementThreshold} outlierQuantile={sampleConfig?.outlierQuantileOverride ?? bundle.config.parameters.outlierQuantile} /></article><article><header><h3>Artefact-cutoff decision</h3><p>Deterministic jitter below the non-outlier quantile with threshold guides.</p></header><ArtefactDecisionPlot families={families} artefactCutoff={summary?.artefactCutoff ?? 0} agreementThreshold={sampleConfig?.agreementOverride ?? bundle.config.parameters.agreementThreshold} outlierQuantile={sampleConfig?.outlierQuantileOverride ?? bundle.config.parameters.outlierQuantile} artefactFraction={sampleConfig?.artefactFractionOverride ?? bundle.config.parameters.artefactFraction} /></article><article><header><h3>Low-agreement positions</h3><p>One longest-run minimum site per non-artefact family, sized by homopolymer run and colored by modal base.</p></header><AgreementPositionPlot consensuses={consensuses} threshold={sampleConfig?.agreementOverride ?? bundle.config.parameters.agreementThreshold} minimumFamilySize={summary?.artefactCutoff ?? 0} /></article><article><header><h3>P(APOBEC|mutations) · classical MDS</h3><p>Pairwise non-gap Hamming distances; point area follows family size and color follows the stored APOBEC posterior.</p></header><MdsApobecPlot records={sampleRecords} /></article><article><header><h3>UMI dinucleotide frequencies</h3><p>Family-unweighted frequencies beside frequencies weighted by reads within each UMI family.</p></header><DinucleotideHeatmaps families={families} /></article></div>
      <TimingSummary bundle={bundle} /><article className="provenance-card"><h3>Audit trail</h3><dl><div><dt>Input SHA-256</dt><dd>{bundle.provenance.inputSha256}</dd></div><div><dt>Config SHA-256</dt><dd>{bundle.provenance.configSha256}</dd></div><div><dt>Engine</dt><dd>{bundle.provenance.engine}</dd></div><div><dt>PORPID source</dt><dd>{bundle.provenance.upstreamBranch}@{bundle.provenance.upstreamCommit.slice(0, 12)}</dd></div>{bundle.inputMappings?.map((mapping) => <div key={`${mapping.slot}-${mapping.uploadedName}`}><dt>{mapping.role} · {mapping.slot}</dt><dd>{mapping.uploadedName}{mapping.expectedName && mapping.expectedName !== mapping.uploadedName ? ` → ${mapping.expectedName}` : ""}</dd></div>)}</dl></article>
    </section>}
    {tab === "families" && <section className="result-section"><div className="table-heading"><div><h2>UMI family decisions</h2><p>Probabilistic offspring assignment, heteroduplex check, length and family-size gates.</p></div><input placeholder="Search UMI, parent, or class" value={familyQuery} onChange={(event) => setFamilyQuery(event.target.value)} /></div><div className="table-scroll"><table><thead><tr><SortHeader label="UMI" column="umi" state={familySort} onChange={setFamilySort} /><SortHeader label="fs" column="familySize" state={familySort} onChange={setFamilySort} /><SortHeader label="classification" column="classification" state={familySort} onChange={setFamilySort} /><SortHeader label="parent" column="mostLikelyParent" state={familySort} onChange={setFamilySort} /><SortHeader label="posterior" column="posteriorProbability" state={familySort} onChange={setFamilySort} /><SortHeader label="minag" column="minimumAgreement" state={familySort} onChange={setFamilySort} /></tr></thead><tbody>{page(filteredFamilies, familyPage).map((row) => <tr key={`${row.sampleIndex}-${row.umi}`}><td><code>{row.umi}</code></td><td>{row.familySize}</td><td><span className={`status ${row.disposition === "likely_real" ? "pass" : "reject"}`}>{row.disposition}</span></td><td><code>{row.mostLikelyParent}</code></td><td>{row.posteriorProbability.toFixed(6)}</td><td>{row.minimumAgreement?.toFixed(2) ?? "—"}</td></tr>)}</tbody></table></div><Pager page={familyPage} count={filteredFamilies.length} onChange={setFamilyPage} /></section>}
    {tab === "sequences" && <section className="result-section"><div className="table-heading"><div><h2>Consensus and post-processing</h2><p>Every filter decision remains inspectable; FASTA exports use the exact stored sequences.</p></div><input placeholder="Search ID or UMI" value={query} onChange={(event) => setQuery(event.target.value)} /></div><div className="table-scroll"><table><thead><tr><SortHeader label="Sequence" column="id" state={sequenceSort} onChange={setSequenceSort} /><SortHeader label="fs" column="familySize" state={sequenceSort} onChange={setSequenceSort} /><SortHeader label="minag" column="minimumAgreement" state={sequenceSort} onChange={setSequenceSort} /><SortHeader label="panel score" column="panelScore" state={sequenceSort} onChange={setSequenceSort} /><SortHeader label="filters" column="filters" state={sequenceSort} onChange={setSequenceSort} /><SortHeader label="APOBEC p" column="apobec" state={sequenceSort} onChange={setSequenceSort} /><SortHeader label="reason" column="reason" state={sequenceSort} onChange={setSequenceSort} /></tr></thead><tbody>{page(records, sequencePage).map((row) => <tr key={row.id}><td><code title={row.id}>{row.id}</code></td><td>{row.familySize}</td><td>{row.minimumAgreement.toFixed(2)}</td><td>{row.panelScore.toFixed(2)}</td><td><span className={`status ${row.artefactPass && row.agreementPass && row.contaminationPass && row.panelPass ? "pass" : "reject"}`}>{row.artefactPass && row.agreementPass && row.contaminationPass && row.panelPass ? "pass" : "reject"}</span></td><td>{row.apobec?.posteriorGaInflated.toFixed(3) ?? "—"}</td><td>{row.rejectionReasons.join("; ") || "—"}</td></tr>)}</tbody></table></div><Pager page={sequencePage} count={records.length} onChange={setSequencePage} /></section>}
    {tab === "contamination" && <section className="result-section"><div className="table-heading"><div><h2>Contamination report</h2><p>Primary calls and the wider zero-proportion suspect pass are retained separately.</p></div><span>{contaminationRows.length.toLocaleString()} calls</span></div><div className="table-scroll"><table><thead><tr><SortHeader label="Sequence" column="sequenceId" state={contaminationSort} onChange={setContaminationSort} /><SortHeader label="Nearest non-self" column="nearestNonselfVariant" state={contaminationSort} onChange={setContaminationSort} /><SortHeader label="distance" column="nearestNonselfDistance" state={contaminationSort} onChange={setContaminationSort} /><SortHeader label="decision" column="decision" state={contaminationSort} onChange={setContaminationSort} /></tr></thead><tbody>{page(contaminationRows, contaminationPage).map((row, index) => <tr key={`${row.sequenceId}-${row.suspectOnly}-${contaminationPage}-${index}`}><td><code>{row.sequenceId}</code></td><td>{row.nearestNonselfVariant}</td><td>{row.nearestNonselfDistance.toPrecision(5)}</td><td><span className={`status ${row.discarded ? "reject" : "warn"}`}>{row.discarded ? "discarded" : row.suspectOnly ? "suspect pass" : "reported"}</span></td></tr>)}</tbody></table></div><Pager page={contaminationPage} count={contaminationRows.length} onChange={setContaminationPage} /></section>}
    {tab === "alignment" && <section className="result-section"><div className="alignment-switch"><div><h2>Tree + alignment workbench</h2><p>Collapsed haplotypes are the default phylogeny. Protein is translated directly from each active nucleotide alignment.</p></div><div><button className={alphabet === "nt" ? "active" : ""} onClick={() => setAlphabet("nt")}>Nucleotide</button><button className={alphabet === "aa" ? "active" : ""} onClick={() => setAlphabet("aa")}>Protein</button></div></div>
      {edited && <div className="alignment-edited-notice"><strong>Edited alignments are present</strong><span>They are stored as separate session copies; original pipeline alignments remain available.</span></div>}{alignmentStatus && <div className="viewer-status">{alignmentStatus}</div>}{alignmentError && <div className="error-box" role="alert">{alignmentError}</div>}
      {alignmentBlock("collapsed")}
      {!showUncollapsed ? <button type="button" className="show-uncollapsed" onClick={() => setShowUncollapsed(true)}>Open uncollapsed family-level phylogeny + alignment</button> : <>{alignmentBlock("uncollapsed")}<button type="button" onClick={() => setShowUncollapsed(false)}>Close uncollapsed view</button></>}
    </section>}
    {tab === "log" && <section className="result-section"><div className="table-heading"><div><h2>Run log</h2><p>Persistent stage summaries, input-slot mappings, interactive tree runs, and fallbacks stored inside the results file.</p></div></div><pre className="run-log">{bundle.log.join("\n")}</pre></section>}
  </main>;
}
