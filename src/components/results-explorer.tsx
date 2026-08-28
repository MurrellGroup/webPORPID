import { useEffect, useMemo, useRef, useState } from "react";
import {
  ALIVIBE_BRIDGE_VERSION, ALIVIBE_SOURCE_REVISION, assertAlivibeInitialLoad, assertAlivibeRoundTripTarget,
  getAlivibeBridge, loadAlivibeNucleotideFasta, readAlivibeNucleotideFasta,
  type AlivibeEditorWindow, type AlivibeMsaJob,
} from "../alivibe-roundtrip";
import { runAlivibeMsa } from "../alivibe-msa-runtime";
import {
  alignmentKey, effectiveAlignment, inspectAlignment, translateAlignmentFasta, validateCorrectedAlignment,
  type AlignmentFrameOffset,
} from "../alignment-utils";
import { runFastTree } from "../biowasm";
import { exportComponent, type ExportKind, safeDatasetName } from "../result-file";
import type { ResultBundle } from "../types";
import { AgreementPositionPlot, ArtefactDecisionPlot, DinucleotideHeatmaps, MdsApobecPlot, UmiDecisionPlot } from "./charts";
import { AlignmentTreeViewer } from "./alignment-tree-viewer";

type Tab = "overview" | "families" | "sequences" | "contamination" | "alignment" | "log";
const PAGE_SIZE = 250;

interface AlivibeSession {
  token: string;
  popup: AlivibeEditorWindow;
  sample: string;
  baseline: string;
  baselineFingerprint: string;
  frameOffset: AlignmentFrameOffset;
  timer?: number;
}

function downloadText(name: string, text: string, mime: string) {
  const url = URL.createObjectURL(new Blob([text], { type: mime })), anchor = document.createElement("a");
  anchor.href = url; anchor.download = name; anchor.click(); setTimeout(() => URL.revokeObjectURL(url), 0);
}

function Pager({ page, count, onChange }: { page: number; count: number; onChange(page: number): void }) {
  const pages = Math.max(1, Math.ceil(count / PAGE_SIZE)), current = Math.min(page, pages - 1);
  if (count <= PAGE_SIZE) return null;
  return <div className="table-pager"><span>Rows {(current * PAGE_SIZE + 1).toLocaleString()}–{Math.min(count, (current + 1) * PAGE_SIZE).toLocaleString()} of {count.toLocaleString()}</span><div><button type="button" disabled={current === 0} onClick={() => onChange(current - 1)}>Previous</button><strong>{current + 1} / {pages}</strong><button type="button" disabled={current + 1 >= pages} onClick={() => onChange(current + 1)}>Next</button></div></div>;
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
  const [kind, setKind] = useState<ExportKind>("consensus-fasta");
  const labels: Array<[ExportKind, string]> = [["consensus-fasta", "Consensus FASTA"], ["passed-consensus-fasta", "Passed consensus FASTA"],
    ["rejected-consensus-fasta", "Rejected consensus FASTA"], ["trimmed-nt-fasta", "Trimmed nucleotide FASTA"], ["trimmed-aa-fasta", "Trimmed amino-acid FASTA"],
    ["family-csv", "UMI family CSV"], ["low-agreement-csv", "Low-agreement CSV"], ["contamination-csv", "Contamination CSV"],
    ["postproc-csv", "Postproc CSV"], ["apobec-csv", "APOBEC CSV"], ["nucleotide-alignment", "Nucleotide alignment"],
    ["protein-alignment", "Translated protein alignment"], ["newick", "Newick tree"], ["log", "Run log"]];
  return <div className="export-menu"><select value={kind} onChange={(event) => setKind(event.target.value as ExportKind)}>{labels.map(([value, label]) => <option value={value} key={value}>{label}</option>)}</select>
    <button type="button" onClick={() => { const result = exportComponent(bundle, kind, sample); downloadText(`${safeDatasetName(bundle.config.dataset)}-${safeDatasetName(sample)}.${result.extension}`, result.text, result.mime); }}>Export</button></div>;
}

function TimingSummary({ bundle }: { bundle: ResultBundle }) {
  if (!bundle.timings?.length) return null;
  return <article className="timing-card"><header><h3>Run timing</h3><p>Wall time measured around each stored pipeline stage.</p></header><div>{bundle.timings.map((entry) => <span key={entry.stage}><strong>{entry.stage.replaceAll("-", " ")}</strong><em>{entry.seconds < 1 ? `${(entry.seconds * 1000).toFixed(1)} ms` : `${entry.seconds.toFixed(2)} s`}</em></span>)}</div></article>;
}

export function ResultsExplorer({ bundle, onSaveResults, onBundleChange }: {
  bundle: ResultBundle; onSaveResults(): void; onBundleChange(bundle: ResultBundle): void;
}) {
  const [tab, setTab] = useState<Tab>("overview"), [sample, setSample] = useState(bundle.summaries[0]?.sample ?? bundle.config.samples[0]?.name ?? "");
  const [query, setQuery] = useState(""), [familyQuery, setFamilyQuery] = useState(""), [alphabet, setAlphabet] = useState<"nt" | "aa">("nt");
  const [familyPage, setFamilyPage] = useState(0), [sequencePage, setSequencePage] = useState(0), [contaminationPage, setContaminationPage] = useState(0);
  const [alignmentStatus, setAlignmentStatus] = useState(""), [alignmentError, setAlignmentError] = useState(""), [alignmentBusy, setAlignmentBusy] = useState(false);
  const bundleRef = useRef(bundle), sampleRef = useRef(sample), changeRef = useRef(onBundleChange), alivibeRef = useRef<AlivibeSession | null>(null);
  bundleRef.current = bundle; sampleRef.current = sample; changeRef.current = onBundleChange;

  const summary = bundle.summaries.find((row) => row.sample === sample);
  const sampleConfig = bundle.config.samples.find((row) => row.name === sample);
  const families = useMemo(() => bundle.umiFamilies.filter((row) => row.sample === sample), [bundle, sample]);
  const consensuses = useMemo(() => bundle.consensuses.filter((row) => row.sample === sample), [bundle, sample]);
  const filteredFamilies = useMemo(() => families.filter((row) => !familyQuery || row.umi.includes(familyQuery.toUpperCase()) || row.mostLikelyParent.includes(familyQuery.toUpperCase()) || row.disposition.toLowerCase().includes(familyQuery.toLowerCase())), [families, familyQuery]);
  const sampleRecords = useMemo(() => bundle.records.filter((row) => row.sample === sample), [bundle, sample]);
  const records = useMemo(() => sampleRecords.filter((row) => !query || row.id.toLowerCase().includes(query.toLowerCase()) || row.umi.includes(query.toUpperCase())), [sampleRecords, query]);
  const contaminationRows = useMemo(() => bundle.contamination.filter((row) => row.sample === sample), [bundle, sample]);
  const effective = useMemo(() => effectiveAlignment(bundle, sample), [bundle, sample]);
  const displayedAlignment = useMemo(() => !effective.fasta ? undefined : alphabet === "nt" ? effective.fasta : translateAlignmentFasta(effective.fasta, effective.frameOffset), [effective, alphabet]);
  const tree = effective.edit ? effective.edit.treeNewick : bundle.trees[alignmentKey(sample)];

  useEffect(() => { setFamilyPage(0); setSequencePage(0); setContaminationPage(0); setAlignmentStatus(""); setAlignmentError(""); }, [sample]);
  useEffect(() => setFamilyPage(0), [familyQuery]); useEffect(() => setSequencePage(0), [query]);
  useEffect(() => () => {
    const session = alivibeRef.current;
    if (session?.timer) window.clearInterval(session.timer);
    if (session && !session.popup.closed) session.popup.close();
  }, []);

  const page = <T,>(rows: T[], index: number) => rows.slice(index * PAGE_SIZE, (index + 1) * PAGE_SIZE);

  async function installCorrection(targetSample: string, baseline: string, corrected: string, frameOffset: AlignmentFrameOffset, source: string) {
    setAlignmentBusy(true); setAlignmentError("");
    try {
      const inspected = validateCorrectedAlignment(baseline, corrected), currentBundle = bundleRef.current, key = alignmentKey(targetSample);
      const original = currentBundle.alignments[key];
      if (!original) throw new Error("The original nucleotide alignment is unavailable.");
      let treeNewick: string | undefined, treeWarning = "";
      try { treeNewick = await runFastTree(inspected.fasta); }
      catch (cause) { treeWarning = ` FastTree could not be refreshed: ${cause instanceof Error ? cause.message : String(cause)}`; }
      const next: ResultBundle = { ...currentBundle, alignmentEdits: { ...currentBundle.alignmentEdits, [key]: {
        fasta: inspected.fasta, frameOffset, baselineFingerprint: inspectAlignment(original).fingerprint,
        editedFingerprint: inspected.fingerprint, source, savedUtc: new Date().toISOString(), treeNewick,
      } } };
      bundleRef.current = next; changeRef.current(next);
      setAlignmentStatus(`Saved corrected alignment: ${inspected.rows.toLocaleString()} rows × ${inspected.columns.toLocaleString()} columns${inspected.removedRows.length ? ` · ${inspected.removedRows.length.toLocaleString()} rows removed` : ""}${inspected.removedNucleotides ? ` · ${inspected.removedNucleotides.toLocaleString()} nucleotide characters removed` : ""} · frame ${frameOffset + 1}.${treeWarning}`);
    } catch (cause) { setAlignmentError(cause instanceof Error ? cause.message : String(cause)); }
    finally { setAlignmentBusy(false); }
  }

  function setFrameOffset(frameOffset: AlignmentFrameOffset) {
    const currentBundle = bundleRef.current, targetSample = sampleRef.current, key = alignmentKey(targetSample), current = effectiveAlignment(currentBundle, targetSample);
    if (!current.fasta) return;
    const fingerprint = inspectAlignment(current.fasta).fingerprint, original = currentBundle.alignments[key] ?? current.fasta;
    const next: ResultBundle = { ...currentBundle, alignmentEdits: { ...currentBundle.alignmentEdits, [key]: {
      fasta: current.fasta, frameOffset, baselineFingerprint: inspectAlignment(original).fingerprint, editedFingerprint: fingerprint,
      source: current.edit?.source ?? "Translation frame selection", savedUtc: new Date().toISOString(), treeNewick: current.edit?.treeNewick ?? currentBundle.trees[key],
    } } };
    bundleRef.current = next; changeRef.current(next); setAlignmentStatus(`Protein translation now starts at nucleotide column ${frameOffset + 1}; this choice is stored in the results file.`);
  }

  function resetCorrection() {
    const currentBundle = bundleRef.current, key = alignmentKey(sampleRef.current), edits = { ...currentBundle.alignmentEdits };
    delete edits[key];
    const next: ResultBundle = { ...currentBundle, alignmentEdits: Object.keys(edits).length ? edits : undefined };
    bundleRef.current = next; changeRef.current(next); setAlignmentError(""); setAlignmentStatus("Restored the pipeline-generated alignment and tree.");
  }

  function createMsaJob(sequences: string[], scoringMode: "nucleotide" | "amino-acid"): AlivibeMsaJob {
    const controller = new AbortController();
    return { result: runAlivibeMsa(sequences, controller.signal, 3, scoringMode), cancel: () => controller.abort() };
  }

  async function returnFromAlivibe(session: AlivibeSession) {
    if (alivibeRef.current?.token !== session.token) throw new Error("This Alivibe editor belongs to an expired round trip.");
    if (session.popup.closed) throw new Error("The Alivibe editor has been closed.");
    const current = effectiveAlignment(bundleRef.current, session.sample).fasta;
    if (!current) throw new Error("The originating nucleotide alignment is no longer loaded.");
    assertAlivibeRoundTripTarget(
      { groupKey: session.sample, alignmentFingerprint: session.baselineFingerprint },
      { groupKey: sampleRef.current, alignmentFingerprint: inspectAlignment(current).fingerprint },
    );
    const returned = readAlivibeNucleotideFasta(session.popup);
    await installCorrection(session.sample, session.baseline, returned.fasta, returned.frameOffset, `Alivibe ${returned.sourceRevision.slice(0, 12)}`);
    alivibeRef.current = null; window.focus(); session.popup.close();
  }

  function openAlivibe() {
    const existing = alivibeRef.current;
    if (existing && !existing.popup.closed) { existing.popup.focus(); setAlignmentError("An Alivibe round trip is already open."); return; }
    const targetSample = sampleRef.current, current = effectiveAlignment(bundleRef.current, targetSample);
    if (!current.fasta) return;
    setAlignmentError(""); setAlignmentStatus("Opening the bundled Alivibe editor…");
    const base = import.meta.env.BASE_URL.endsWith("/") ? import.meta.env.BASE_URL : `${import.meta.env.BASE_URL}/`;
    const editorUrl = new URL(`${base}tools/alivibe.html`, window.location.origin);
    editorUrl.searchParams.set("swigBridge", String(ALIVIBE_BRIDGE_VERSION)); editorUrl.searchParams.set("source", ALIVIBE_SOURCE_REVISION.slice(0, 12));
    const token = globalThis.crypto?.randomUUID?.() ?? `${Date.now()}-${Math.random().toString(36).slice(2)}`;
    const popup = window.open(editorUrl.href, `webporpid-alivibe-${token}`, "popup,width=1500,height=920") as AlivibeEditorWindow | null;
    if (!popup) { setAlignmentError("The browser blocked the Alivibe window. Allow pop-ups for this site and try again."); return; }
    const session: AlivibeSession = { token, popup, sample: targetSample, baseline: current.fasta, baselineFingerprint: inspectAlignment(current.fasta).fingerprint, frameOffset: current.frameOffset };
    alivibeRef.current = session; let attempts = 0;
    session.timer = window.setInterval(() => {
      attempts += 1;
      if (popup.closed) { window.clearInterval(session.timer); if (alivibeRef.current?.token === token) alivibeRef.current = null; setAlignmentStatus("Alivibe closed without returning an alignment."); return; }
      try {
        const bridge = getAlivibeBridge(popup);
        if (!bridge) { if (attempts < 240) return; throw new Error("The bundled Alivibe bridge did not become ready."); }
        bridge.installMsaRunner(createMsaJob);
        bridge.installFastTreeRunner(runFastTree);
        const loaded = loadAlivibeNucleotideFasta(popup, session.baseline, session.frameOffset); assertAlivibeInitialLoad(session.baseline, loaded);
        const controls = popup.document.getElementById("controls"); if (!controls) return;
        if (!popup.document.getElementById("webporpid-return-control")) {
          const group = popup.document.createElement("div"); group.id = "webporpid-return-control"; group.className = "control-group";
          const label = popup.document.createElement("label"); label.textContent = "webPORPID round trip";
          const button = popup.document.createElement("button"); button.type = "button"; button.textContent = "Return alignment to webPORPID"; button.className = "active";
          button.onclick = () => void returnFromAlivibe(session).catch((cause) => { setAlignmentError(cause instanceof Error ? cause.message : String(cause)); window.focus(); });
          group.append(label, button); controls.prepend(group);
        }
        setAlignmentStatus(`Exact nucleotide alignment loaded in bundled Alivibe ${ALIVIBE_SOURCE_REVISION.slice(0, 12)}. Edit it, then press Return alignment to webPORPID.`);
        window.clearInterval(session.timer);
      } catch (cause) {
        setAlignmentError(cause instanceof Error ? cause.message : String(cause)); setAlignmentStatus("Alivibe round trip stopped before editing.");
        if (alivibeRef.current?.token === token) alivibeRef.current = null; window.clearInterval(session.timer);
      }
    }, 250);
  }

  async function importCorrected(file?: File) {
    if (!file || !effective.fasta) return;
    await installCorrection(sample, effective.fasta, await file.text(), effective.frameOffset, `Imported corrected FASTA · ${file.name}`);
  }

  return <main className="results-page">
    <section className="results-hero"><div><span className="section-kicker">Loaded analysis</span><h1>{bundle.config.dataset}</h1><p>{bundle.provenance.inputName} · {bundle.provenance.createdUtc} · {bundle.provenance.workers} workers</p></div>
      <div className="result-actions"><select aria-label="Sample" value={sample} onChange={(event) => setSample(event.target.value)}>{bundle.summaries.map((row) => <option key={row.sample}>{row.sample}</option>)}</select><ExportMenu bundle={bundle} sample={sample} /><button className="primary" type="button" onClick={onSaveResults}>Save results file</button></div></section>
    <nav className="result-tabs">{(["overview", "families", "sequences", "contamination", "alignment", "log"] as Tab[]).map((value) => <button type="button" className={tab === value ? "active" : ""} onClick={() => setTab(value)} key={value}>{value}</button>)}</nav>
    {tab === "overview" && <section className="result-section">
      <div className="metric-grid"><article><span>Demultiplexed</span><strong>{summary?.demultiplexedReads.toLocaleString() ?? "0"}</strong></article><article><span>Observed UMIs</span><strong>{summary?.observedUmis.toLocaleString() ?? "0"}</strong></article><article><span>Consensus</span><strong>{summary?.consensusSequences.toLocaleString() ?? "0"}</strong></article><article><span>Postproc passed</span><strong>{summary?.postprocPassed.toLocaleString() ?? "0"}</strong></article><article><span>Artefact cutoff</span><strong>{summary?.artefactCutoff ?? 0}</strong></article></div>
      <Filters bundle={bundle} sample={sample} />
      <div className="chart-grid report-figures">
        <article><header><h3>UMI family size × UMI length</h3><p>Julia report parity: PORPID class, family-size and UMI-length decisions with quantile and artefact thresholds.</p></header><UmiDecisionPlot families={families} artefactCutoff={summary?.artefactCutoff ?? 0} agreementThreshold={sampleConfig?.agreementOverride ?? bundle.config.parameters.agreementThreshold} outlierQuantile={sampleConfig?.outlierQuantileOverride ?? bundle.config.parameters.outlierQuantile} /></article>
        <article><header><h3>Artefact-cutoff decision</h3><p>Julia report parity: deterministic jitter below the non-outlier quantile and the 5%–95% guide thresholds.</p></header><ArtefactDecisionPlot families={families} artefactCutoff={summary?.artefactCutoff ?? 0} agreementThreshold={sampleConfig?.agreementOverride ?? bundle.config.parameters.agreementThreshold} outlierQuantile={sampleConfig?.outlierQuantileOverride ?? bundle.config.parameters.outlierQuantile} artefactFraction={sampleConfig?.artefactFractionOverride ?? bundle.config.parameters.artefactFraction} /></article>
        <article><header><h3>Low-agreement positions</h3><p>Julia report parity: one longest-run minimum site per non-artefact family, sized by homopolymer run and colored by modal base.</p></header><AgreementPositionPlot consensuses={consensuses} threshold={sampleConfig?.agreementOverride ?? bundle.config.parameters.agreementThreshold} minimumFamilySize={summary?.artefactCutoff ?? 0} /></article>
        <article><header><h3>P(APOBEC|mutations) · classical MDS</h3><p>Pairwise non-gap Hamming distances; point area follows family size and color follows the stored APOBEC posterior.</p></header><MdsApobecPlot records={sampleRecords} /></article>
        <article><header><h3>UMI dinucleotide frequencies</h3><p>Julia report parity: family-unweighted frequencies beside frequencies weighted by the number of reads in each UMI family.</p></header><DinucleotideHeatmaps families={families} /></article>
      </div>
      <TimingSummary bundle={bundle} />
      <article className="provenance-card"><h3>Audit trail</h3><dl><div><dt>Input SHA-256</dt><dd>{bundle.provenance.inputSha256}</dd></div><div><dt>Config SHA-256</dt><dd>{bundle.provenance.configSha256}</dd></div><div><dt>Engine</dt><dd>{bundle.provenance.engine}</dd></div><div><dt>PORPID source</dt><dd>{bundle.provenance.upstreamBranch}@{bundle.provenance.upstreamCommit.slice(0, 12)}</dd></div></dl></article>
    </section>}
    {tab === "families" && <section className="result-section"><div className="table-heading"><div><h2>UMI family decisions</h2><p>Probabilistic offspring assignment, heteroduplex check, length and family-size gates.</p></div><input placeholder="Search UMI, parent, or class" value={familyQuery} onChange={(event) => setFamilyQuery(event.target.value)} /></div><div className="table-scroll"><table><thead><tr><th>UMI</th><th>fs</th><th>classification</th><th>parent</th><th>posterior</th><th>minag</th></tr></thead><tbody>{page(filteredFamilies, familyPage).map((row) => <tr key={`${row.sampleIndex}-${row.umi}`}><td><code>{row.umi}</code></td><td>{row.familySize}</td><td><span className={`status ${row.disposition === "likely_real" ? "pass" : "reject"}`}>{row.disposition}</span></td><td><code>{row.mostLikelyParent}</code></td><td>{row.posteriorProbability.toFixed(6)}</td><td>{row.minimumAgreement?.toFixed(2) ?? "—"}</td></tr>)}</tbody></table></div><Pager page={familyPage} count={filteredFamilies.length} onChange={setFamilyPage} /></section>}
    {tab === "sequences" && <section className="result-section"><div className="table-heading"><div><h2>Consensus and post-processing</h2><p>Every filter decision remains inspectable; FASTA exports use the exact stored sequences.</p></div><input placeholder="Search ID or UMI" value={query} onChange={(event) => setQuery(event.target.value)} /></div><div className="table-scroll"><table><thead><tr><th>Sequence</th><th>fs</th><th>minag</th><th>panel score</th><th>filters</th><th>APOBEC p</th><th>reason</th></tr></thead><tbody>{page(records, sequencePage).map((row) => <tr key={row.id}><td><code title={row.id}>{row.id}</code></td><td>{row.familySize}</td><td>{row.minimumAgreement.toFixed(2)}</td><td>{row.panelScore.toFixed(2)}</td><td><span className={`status ${row.artefactPass && row.agreementPass && row.contaminationPass && row.panelPass ? "pass" : "reject"}`}>{row.artefactPass && row.agreementPass && row.contaminationPass && row.panelPass ? "pass" : "reject"}</span></td><td>{row.apobec?.posteriorGaInflated.toFixed(3) ?? "—"}</td><td>{row.rejectionReasons.join("; ") || "—"}</td></tr>)}</tbody></table></div><Pager page={sequencePage} count={records.length} onChange={setSequencePage} /></section>}
    {tab === "contamination" && <section className="result-section"><div className="table-heading"><div><h2>Contamination report</h2><p>Primary calls and the wider zero-proportion suspect pass are retained separately.</p></div><span>{contaminationRows.length.toLocaleString()} calls</span></div><div className="table-scroll"><table><thead><tr><th>Sequence</th><th>Nearest non-self</th><th>distance</th><th>decision</th></tr></thead><tbody>{page(contaminationRows, contaminationPage).map((row, index) => <tr key={`${row.sequenceId}-${row.suspectOnly}-${contaminationPage}-${index}`}><td><code>{row.sequenceId}</code></td><td>{row.nearestNonselfVariant}</td><td>{row.nearestNonselfDistance.toPrecision(5)}</td><td><span className={`status ${row.discarded ? "reject" : "warn"}`}>{row.discarded ? "discarded" : row.suspectOnly ? "suspect pass" : "reported"}</span></td></tr>)}</tbody></table></div><Pager page={contaminationPage} count={contaminationRows.length} onChange={setContaminationPage} /></section>}
    {tab === "alignment" && <section className="result-section">
      <div className="alignment-switch"><div><h2>Tree + alignment workbench</h2><p>Swig’s branch-length layout is coordinated with a virtualized alignment; protein is translated directly from the active nucleotide rows.</p></div><div><button className={alphabet === "nt" ? "active" : ""} onClick={() => setAlphabet("nt")}>Nucleotide</button><button className={alphabet === "aa" ? "active" : ""} onClick={() => setAlphabet("aa")}>Protein</button></div></div>
      {effective.fasta && <div className="alignment-edit-bar"><button type="button" disabled={alignmentBusy} onClick={openAlivibe}>Open in Alivibe ↗</button><label className="button-like">Import corrected FASTA<input type="file" accept=".fasta,.fa,.fas,.fna,text/plain" onChange={(event) => void importCorrected(event.target.files?.[0])} /></label><label><span>Protein frame</span><select value={effective.frameOffset} onChange={(event) => setFrameOffset(Number(event.target.value) as AlignmentFrameOffset)}><option value="0">Start at nucleotide column 1</option><option value="1">Start at nucleotide column 2</option><option value="2">Start at nucleotide column 3</option></select></label>{effective.edit && <button type="button" onClick={resetCorrection}>Restore pipeline alignment</button>}<span>{effective.edit ? `Edited · ${effective.edit.source}` : "Pipeline alignment"}</span></div>}
      {alignmentStatus && <div className="viewer-status">{alignmentStatus}</div>}{alignmentError && <div className="error-box" role="alert">{alignmentError}</div>}
      {displayedAlignment ? <AlignmentTreeViewer fasta={displayedAlignment} newick={tree} alphabet={alphabet} name={`${safeDatasetName(bundle.config.dataset)}-${safeDatasetName(sample)}`} /> : <div className="empty-state">No nucleotide alignment passed the selected sample’s non-functional filters.</div>}
    </section>}
    {tab === "log" && <section className="result-section"><div className="table-heading"><div><h2>Run log</h2><p>Persistent stage summaries and fallbacks stored inside the results file.</p></div></div><pre className="run-log">{bundle.log.join("\n")}</pre></section>}
  </main>;
}
