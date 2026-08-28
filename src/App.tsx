import { useMemo, useRef, useState, type DragEvent } from "react";
import { blankConfig, parseConfigYaml, resolveReferenceFiles, serializeConfigYaml } from "./config";
import { nameMatchingSlot, referenceFileMap, referenceMappingRecords, referenceSlots, type ReferenceSlot } from "./input-mapping";
import { decodeResult, encodeResult, safeDatasetName } from "./result-file";
import type { InputFileMapping, PipelineConfig, PipelineProgress, ResultBundle } from "./types";
import { ResultsExplorer } from "./components/results-explorer";
import { ConfigForm } from "./components/config-form";

function download(name: string, data: string | Uint8Array, mime: string) {
  const body: BlobPart = typeof data === "string" ? data : data.slice().buffer as ArrayBuffer;
  const url = URL.createObjectURL(new Blob([body], { type: mime }));
  const anchor = document.createElement("a");
  anchor.href = url; anchor.download = name; anchor.click();
  setTimeout(() => URL.revokeObjectURL(url), 0);
}

interface SlotConflict { file: File; target: ReferenceSlot; matching: ReferenceSlot }

const fileKey = (file: File) => `${file.name}\0${file.size}\0${file.lastModified}`;
const isYaml = (file: File) => /\.ya?ml$/i.test(file.name);
const isFastq = (file: File) => /\.(?:fastq|fq)(?:\.gz)?$/i.test(file.name);
const isResult = (file: File) => /\.webporpid$/i.test(file.name);

function StageProgress({ value }: { value: PipelineProgress }) {
  const percent = Math.max(0, Math.min(100, Math.round(value.fraction * 100)));
  return <section className="run-progress" aria-live="polite">
    <div><span>{value.stage}</span><strong>{percent}%</strong></div>
    <progress max="100" value={percent} />
    <p>{value.detail}</p>
  </section>;
}

export default function App() {
  const defaultWorkers = Math.max(1, navigator.hardwareConcurrency || 1);
  const [draftConfig, setDraftConfig] = useState<PipelineConfig>(() => blankConfig());
  const [fastq, setFastq] = useState<File>();
  const [configFile, setConfigFile] = useState<File>();
  const [configName, setConfigName] = useState("new configuration");
  const [configText, setConfigText] = useState(() => serializeConfigYaml(blankConfig()));
  const [configIssue, setConfigIssue] = useState("");
  const [referenceAssignments, setReferenceAssignments] = useState<Record<string, File>>({});
  const [unassignedReferences, setUnassignedReferences] = useState<File[]>([]);
  const [slotConflict, setSlotConflict] = useState<SlotConflict>();
  const [workers, setWorkers] = useState(defaultWorkers);
  const [deferPhylogeny, setDeferPhylogeny] = useState(false);
  const [progress, setProgress] = useState<PipelineProgress>();
  const [result, setResult] = useState<ResultBundle>();
  const [error, setError] = useState("");
  const [running, setRunning] = useState(false);
  const [cancelling, setCancelling] = useState(false);
  const workerRef = useRef<Worker | undefined>(undefined);
  const maxWorkers = useMemo(() => Math.max(1, navigator.hardwareConcurrency || 1), []);
  const slots = useMemo(() => referenceSlots(draftConfig), [draftConfig]);

  function reconcileReferences(config: PipelineConfig, extras: File[] = []) {
    const nextSlots = referenceSlots(config), valid = new Set(nextSlots.map((slot) => slot.id));
    const next: Record<string, File> = {};
    const pool: File[] = [...unassignedReferences, ...extras];
    for (const [id, file] of Object.entries(referenceAssignments)) {
      if (valid.has(id)) next[id] = file;
      else pool.push(file);
    }
    const leftovers: File[] = [];
    for (const file of pool) {
      const match = nameMatchingSlot(file.name, nextSlots);
      if (match && !next[match.id]) next[match.id] = file;
      else if (!Object.values(next).some((candidate) => fileKey(candidate) === fileKey(file)) && !leftovers.some((candidate) => fileKey(candidate) === fileKey(file))) leftovers.push(file);
    }
    setReferenceAssignments(next); setUnassignedReferences(leftovers);
  }

  async function chooseConfig(file?: File, referenceFiles: File[] = []): Promise<PipelineConfig | undefined> {
    if (!file) return undefined;
    const source = await file.text();
    setConfigFile(file); setConfigName(file.name); setConfigText(source); setError("");
    try {
      const parsed = parseConfigYaml(source); setDraftConfig(parsed); setConfigIssue(""); reconcileReferences(parsed, referenceFiles); return parsed;
    } catch (cause) { setConfigIssue(cause instanceof Error ? cause.message : String(cause)); return undefined; }
  }

  function changeConfig(config: PipelineConfig) {
    setConfigFile(undefined); setDraftConfig(config); setConfigName("edited configuration"); reconcileReferences(config);
    try { setConfigText(serializeConfigYaml(config)); setConfigIssue(""); }
    catch (cause) { setConfigIssue(cause instanceof Error ? cause.message : String(cause)); }
  }

  function changeRawYaml(source: string) {
    setConfigFile(undefined); setConfigText(source); setConfigName("edited configuration");
    try { const parsed = parseConfigYaml(source); setDraftConfig(parsed); setConfigIssue(""); reconcileReferences(parsed); }
    catch (cause) { setConfigIssue(cause instanceof Error ? cause.message : String(cause)); }
  }

  function startFromScratch() {
    const config = blankConfig();
    setConfigFile(undefined); setDraftConfig(config); setConfigText(serializeConfigYaml(config)); setConfigName("new configuration"); setConfigIssue(""); reconcileReferences(config);
  }

  function addReferenceFiles(files: File[], target?: ReferenceSlot, config = draftConfig) {
    const nextSlots = referenceSlots(config);
    if (target && files[0]) {
      const matching = nameMatchingSlot(files[0].name, nextSlots, target.id);
      if (matching) { setSlotConflict({ file: files[0], target, matching }); return; }
    }
    const next = { ...referenceAssignments }, loose = [...unassignedReferences];
    const retain = (file: File | undefined) => {
      if (file && !Object.values(next).some((candidate) => fileKey(candidate) === fileKey(file)) && !loose.some((candidate) => fileKey(candidate) === fileKey(file))) loose.push(file);
    };
    const assign = (slot: ReferenceSlot, file: File) => {
      const previous = next[slot.id]; next[slot.id] = file;
      const index = loose.findIndex((candidate) => fileKey(candidate) === fileKey(file)); if (index >= 0) loose.splice(index, 1);
      if (previous && fileKey(previous) !== fileKey(file)) retain(previous);
    };
    for (const [index, file] of files.entries()) {
      const slot = index === 0 && target ? target : nameMatchingSlot(file.name, nextSlots);
      if (slot) assign(slot, file);
      else retain(file);
    }
    setReferenceAssignments(next); setUnassignedReferences(loose);
  }

  async function routeFiles(files: File[]) {
    let config = draftConfig;
    const resultFile = files.find(isResult); if (resultFile) { await loadResults(resultFile); return; }
    const reads = files.filter(isFastq); if (reads.length) setFastq(reads.at(-1));
    const references = files.filter((file) => !isYaml(file) && !isFastq(file) && !isResult(file));
    const yaml = files.find(isYaml);
    if (yaml) config = await chooseConfig(yaml, references) ?? config;
    else if (references.length) addReferenceFiles(references, undefined, config);
  }

  function dropFiles(event: DragEvent<HTMLElement>, handler: (files: File[]) => void | Promise<void>) {
    event.preventDefault(); const files = [...event.dataTransfer.files]; if (files.length) void handler(files);
  }

  function resolveSlotConflict(action: "matching" | "chosen" | "swap") {
    if (!slotConflict) return;
    const { file, target, matching } = slotConflict, next = { ...referenceAssignments }, loose = unassignedReferences.filter((entry) => fileKey(entry) !== fileKey(file));
    const retain = (candidate?: File) => {
      if (candidate && fileKey(candidate) !== fileKey(file) && !Object.values(next).some((assigned) => fileKey(assigned) === fileKey(candidate)) && !loose.some((entry) => fileKey(entry) === fileKey(candidate))) loose.push(candidate);
    };
    if (action === "chosen") { const previous = next[target.id]; next[target.id] = file; retain(previous); }
    else if (action === "matching") { const previous = next[matching.id]; next[matching.id] = file; retain(previous); }
    else {
      const previousMatch = next[matching.id], previousTarget = next[target.id]; next[matching.id] = file;
      if (previousMatch) next[target.id] = previousMatch; else delete next[target.id];
      retain(previousTarget);
    }
    setReferenceAssignments(next); setUnassignedReferences(loose); setSlotConflict(undefined);
  }

  async function run() {
    if (!fastq) { setError("Choose a FASTQ or FASTQ.GZ input file."); return; }
    if (!configText.trim()) { setError("Choose or paste a PORPID YAML configuration."); return; }
    setError(""); setResult(undefined); setRunning(true); setCancelling(false);
    try {
      const missing = slots.filter((slot) => !referenceAssignments[slot.id]);
      if (missing.length) throw new Error(`Assign a file to every required reference slot: ${missing.map((slot) => slot.expectedName).join(", ")}.`);
      const config = await resolveReferenceFiles(parseConfigYaml(configText), referenceFileMap(referenceAssignments));
      const inputMappings: InputFileMapping[] = [
        { slot: "reads", role: "reads", uploadedName: fastq.name, uploadedSize: fastq.size },
        { slot: "configuration", role: "configuration", uploadedName: configFile?.name ?? configName, uploadedSize: configFile?.size ?? new TextEncoder().encode(configText).byteLength },
        ...referenceMappingRecords(slots, referenceAssignments),
      ];
      const worker = new Worker(new URL("./pipeline-worker.ts", import.meta.url), { type: "module" });
      workerRef.current = worker;
      worker.onmessage = (event: MessageEvent<{ type: "progress"; progress: PipelineProgress } | { type: "result"; result: ResultBundle } | { type: "error"; message: string }>) => {
        if (event.data.type === "progress") setProgress(event.data.progress);
        if (event.data.type === "result") { setResult(event.data.result); setRunning(false); setCancelling(false); worker.terminate(); workerRef.current = undefined; }
        if (event.data.type === "error") { setError(event.data.message); setRunning(false); setCancelling(false); worker.terminate(); workerRef.current = undefined; }
      };
      worker.onerror = (event) => { setError(event.message || "The pipeline worker failed."); setRunning(false); setCancelling(false); worker.terminate(); workerRef.current = undefined; };
      worker.postMessage({ type: "run", file: fastq, config, workers, deferPhylogeny, inputMappings });
    } catch (cause) { setError(cause instanceof Error ? cause.message : String(cause)); setRunning(false); setCancelling(false); }
  }

  function cancel() {
    workerRef.current?.postMessage({ type: "cancel" });
    setCancelling(true); setError("");
    setProgress((current) => current ? { ...current, detail: "Cancelling safely and removing temporary partitions…" } : current);
  }

  async function loadResults(file?: File) {
    if (!file) return;
    setError("");
    try { setResult(decodeResult(new Uint8Array(await file.arrayBuffer()))); }
    catch (cause) { setError(cause instanceof Error ? cause.message : String(cause)); }
  }

  function saveResults() {
    if (!result) return;
    download(`${safeDatasetName(result.config.dataset)}.webporpid`, encodeResult(result), "application/vnd.webporpid.results");
  }

  function saveConfig() {
    try {
      download(`${safeDatasetName(draftConfig.dataset)}.yaml`, serializeConfigYaml(draftConfig), "application/yaml");
      setConfigIssue("");
    } catch (cause) { setConfigIssue(cause instanceof Error ? cause.message : String(cause)); }
  }

  return <div className="app-shell">
    <header className="site-header"><a className="brand" href={import.meta.env.BASE_URL} aria-label="webPORPID home"><span className="brand-mark">wp</span><span>webPORPID<small>Nanopore &amp; PacBio analysis</small></span></a><nav><a href="#run">Run</a><a href="#about">Methods</a><a href="https://github.com/MurrellGroup/webPORPID">GitHub</a></nav></header>
    {!result ? <>
      <main className="landing">
        <section className="intro"><div><span className="section-kicker">PORPID, entirely on your machine</span><h1>From long reads to auditable within-host variants.</h1><p>Stream compressed FASTQ through demultiplexing, probabilistic UMI grouping, indel-aware consensus, contamination control, panel and functional filters—without uploading sequence data.</p><div className="privacy-pill"><span />Input reads stay in this browser or CLI process.</div></div><div className="intro-card"><strong>C++20 · WASI · SIMD</strong><p>The same portable core runs in browser workers and <code>porpid-cli</code>. Intermediate reads are partitioned to disk-backed storage and released after consensus.</p><dl><div><dt>Input</dt><dd>.fastq / .fastq.gz</dd></div><div><dt>Output</dt><dd>.webporpid</dd></div><div><dt>Default CPUs</dt><dd>{maxWorkers}</dd></div></dl></div></section>
        <section className="runner" id="run"><div className="section-heading"><div><span className="section-kicker">New analysis</span><h2>Configure the run</h2></div><div className="runner-shortcuts"><label className="secondary button-like">Load results<input type="file" accept=".webporpid" onChange={(event) => void loadResults(event.target.files?.[0])} /></label></div></div>
          <label className="all-files-drop" onDragOver={(event) => event.preventDefault()} onDrop={(event) => dropFiles(event, routeFiles)}><strong>Drop any run files here</strong><span>Files accumulate across drops. YAML, FASTQ, reference FASTAs, and saved results are routed automatically.</span><input type="file" multiple onChange={(event) => { const files = [...event.target.files ?? []]; event.currentTarget.value = ""; void routeFiles(files); }} /></label>
          <div className="input-grid">
            <label className="file-card" onDragOver={(event) => event.preventDefault()} onDrop={(event) => dropFiles(event, (files) => { const reads = files.filter(isFastq); if (reads.length) setFastq(reads.at(-1)); })}><span className="step">1</span><strong>Compressed reads</strong><small>FASTQ or FASTQ.GZ; decompressed incrementally.</small><input type="file" accept=".fastq,.fq,.gz,application/gzip" onChange={(event) => { setFastq(event.target.files?.[0]); event.currentTarget.value = ""; }} /><em>{fastq ? `${fastq.name} · ${(fastq.size / 1048576).toFixed(1)} MiB` : "Choose or drop reads"}</em></label>
            <label className="file-card" onDragOver={(event) => event.preventDefault()} onDrop={(event) => dropFiles(event, async (files) => { const yaml = files.find(isYaml); if (yaml) await chooseConfig(yaml); })}><span className="step">2</span><strong>PORPID configuration</strong><small>Original nanopore-branch YAML is accepted.</small><input type="file" accept=".yaml,.yml,text/yaml" onChange={(event) => { const file = event.target.files?.[0]; event.currentTarget.value = ""; void chooseConfig(file); }} /><em>{configName || "Choose or drop YAML"}</em></label>
            <label className="file-card" onDragOver={(event) => event.preventDefault()} onDrop={(event) => dropFiles(event, (files) => addReferenceFiles(files))}><span className="step">3</span><strong>Reference FASTAs</strong><small>Add files in one selection or one after another. Exact filenames are assigned automatically.</small><input type="file" accept=".fasta,.fa,.fas,.fna,text/plain" multiple onChange={(event) => { const files = [...event.target.files ?? []]; event.currentTarget.value = ""; addReferenceFiles(files); }} /><em>{Object.keys(referenceAssignments).length || unassignedReferences.length ? `${Object.keys(referenceAssignments).length} assigned · ${unassignedReferences.length} unassigned` : "Choose or drop references"}</em></label>
          </div>
          <section className="reference-mapper"><header><div><span className="section-kicker">YAML-aware file mapping</span><h3>Required reference slots</h3><p>Each slot is named from the active configuration. A deliberately renamed upload can still be assigned to its intended slot.</p></div><span>{slots.filter((slot) => referenceAssignments[slot.id]).length}/{slots.length} assigned</span></header>
            <div className="reference-slot-grid">{slots.map((slot) => { const assigned = referenceAssignments[slot.id]; return <div className={`reference-slot ${assigned ? "assigned" : ""}`} key={slot.id} onDragOver={(event) => event.preventDefault()} onDrop={(event) => dropFiles(event, (files) => addReferenceFiles(files, slot))}><label><strong>{slot.label}{slot.samples.length ? ` · ${slot.samples.join(", ")}` : ""}</strong><small>Expected: {slot.expectedName}</small><em>{assigned ? `${assigned.name} · ${(assigned.size / 1024).toFixed(1)} KiB` : "Drop a FASTA here"}</em><input type="file" accept=".fasta,.fa,.fas,.fna,text/plain" onChange={(event) => { const file = event.target.files?.[0]; event.currentTarget.value = ""; if (file) addReferenceFiles([file], slot); }} /></label>{assigned && <button type="button" aria-label={`Remove ${assigned.name} from ${slot.label}`} onClick={() => { const next = { ...referenceAssignments }; delete next[slot.id]; setReferenceAssignments(next); setUnassignedReferences((current) => [...current, assigned]); }}>×</button>}</div>; })}</div>
            {unassignedReferences.length > 0 && <div className="unassigned-files"><strong>Unassigned files</strong>{unassignedReferences.map((file) => <div key={fileKey(file)}><span>{file.name}</span><select defaultValue="" onChange={(event) => { const slot = slots.find((candidate) => candidate.id === event.target.value); if (!slot) return; addReferenceFiles([file], slot); setUnassignedReferences((current) => current.filter((candidate) => fileKey(candidate) !== fileKey(file))); }}><option value="">Assign to…</option>{slots.map((slot) => <option value={slot.id} key={slot.id}>{slot.label} · {slot.expectedName}</option>)}</select><button type="button" onClick={() => setUnassignedReferences((current) => current.filter((candidate) => fileKey(candidate) !== fileKey(file)))}>Remove</button></div>)}</div>}
          </section>
          <ConfigForm config={draftConfig} onChange={changeConfig} onReset={startFromScratch} onDownload={saveConfig} />
          <details className="config-editor"><summary><span>Raw YAML</span><small>{configIssue ? "Needs attention" : "Synchronized with the form"}</small></summary><textarea value={configText} onChange={(event) => changeRawYaml(event.target.value)} spellCheck={false} placeholder="dataset:\n  samples:\n    ..." /></details>
          {configIssue && <div className="config-warning" role="status"><strong>Configuration is not ready</strong><span>{configIssue}</span><small>Continue editing in the form or raw YAML. The analysis will not start until the YAML is valid.</small></div>}
          <div className="run-bar"><label><span>CPU workers</span><input type="number" min="1" max={maxWorkers} value={workers} onChange={(event) => setWorkers(Math.max(1, Math.min(maxWorkers, Number(event.target.value) || 1)))} /><small>Maximum detected: {maxWorkers}</small></label><label className="defer-tree"><input type="checkbox" checked={deferPhylogeny} onChange={(event) => setDeferPhylogeny(event.target.checked)} /><span>Defer phylogeny inference</span><small>Store collapsed alignments now and infer trees later from the results explorer.</small></label><div>{running ? <button className="danger" type="button" disabled={cancelling} onClick={cancel}>{cancelling ? "Cancelling…" : "Cancel analysis"}</button> : <button className="primary run-button" type="button" onClick={() => void run()}>Run webPORPID</button>}</div></div>
          {progress && running && <StageProgress value={progress} />}
          {error && <div className="error-box" role="alert">{error}</div>}
        </section>
        <section className="method-strip" id="about"><article><span>01</span><h3>Stream &amp; demultiplex</h3><p>Gzip chunks are decoded incrementally. Read-quality, primer, orientation, sample-ID and BPB logic follows the nanopore branch.</p></article><article><span>02</span><h3>Group &amp; call consensus</h3><p>Sparse two-error offspring likelihoods, LDA decisions, heteroduplex QC, seeded alignment and minimum-agreement counting.</p></article><article><span>03</span><h3>Filter &amp; explore</h3><p>Run-aware contamination, panel and functional filters, APOBEC model, aligned variants, phylogeny and component exports.</p></article></section>
      </main>
    </> : <ResultsExplorer bundle={result} onSaveResults={saveResults} onBundleChange={setResult} />}
    {result && <button className="new-analysis" type="button" onClick={() => { setResult(undefined); setProgress(undefined); }}>← New or load another analysis</button>}
    {slotConflict && <div className="modal-backdrop" role="presentation"><section className="slot-conflict" role="dialog" aria-modal="true" aria-labelledby="slot-conflict-title"><span className="section-kicker">Filename conflict</span><h3 id="slot-conflict-title">This file matches a different YAML slot</h3><p><strong>{slotConflict.file.name}</strong> was dropped onto <strong>{slotConflict.target.expectedName}</strong>, but its filename matches <strong>{slotConflict.matching.expectedName}</strong>.</p><div><button type="button" className="primary" onClick={() => resolveSlotConflict("swap")}>Swap into filename-matching slot</button><button type="button" onClick={() => resolveSlotConflict("matching")}>Use matching slot</button><button type="button" onClick={() => resolveSlotConflict("chosen")}>Keep chosen slot</button><button type="button" onClick={() => setSlotConflict(undefined)}>Cancel</button></div></section></div>}
    <footer className="site-footer"><span>webPORPID · standalone local analysis</span><span>Deterministic, inspectable processing and session files</span></footer>
  </div>;
}
