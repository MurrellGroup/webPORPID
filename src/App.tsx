import { useEffect, useMemo, useRef, useState, type DragEvent } from "react";
import { blankConfig, parseConfigYaml, resolveReferenceFiles, serializeConfigYaml } from "./config";
import { nameMatchingSlot, referenceFileMap, referenceMappingRecords, referenceSlots, type ReferenceSlot } from "./input-mapping";
import type { ExternalScratchDirectoryHandle } from "./partition-store";
import { decodeResult, encodeResult, safeDatasetName } from "./result-file";
import type { InputFileMapping, OptionalStageName, PipelineConfig, PipelineProgress, ResultBundle } from "./types";
import { ResultsExplorer } from "./components/results-explorer";
import { ConfigForm } from "./components/config-form";
import { MethodLink } from "./components/method-link";
import packageInformation from "../package.json";

export const APP_VERSION = packageInformation.version;

function download(name: string, data: string | Uint8Array, mime: string) {
  const body: BlobPart = typeof data === "string" ? data : data.slice().buffer as ArrayBuffer;
  const url = URL.createObjectURL(new Blob([body], { type: mime }));
  const anchor = document.createElement("a");
  anchor.href = url; anchor.download = name; anchor.click();
  setTimeout(() => URL.revokeObjectURL(url), 0);
}

interface SlotConflict { file: File; target: ReferenceSlot; matching: ReferenceSlot }
type SpoolStorage = "automatic" | "external-directory";
type DirectoryPickerWindow = Window & { showDirectoryPicker?: (options?: { id?: string; mode?: "read" | "readwrite" }) => Promise<ExternalScratchDirectoryHandle> };

const fileKey = (file: File) => `${file.name}\0${file.size}\0${file.lastModified}`;
const isYaml = (file: File) => /\.ya?ml$/i.test(file.name);
const isFastq = (file: File) => /\.(?:fastq|fq)(?:\.gz)?$/i.test(file.name);
const isResult = (file: File) => /\.webporpid$/i.test(file.name);

function StageProgress({ value, onSkip, skipping }: { value: PipelineProgress; onSkip?(stage: OptionalStageName): void; skipping?: boolean }) {
  const [clock, setClock] = useState(() => Date.now()), started = useRef(Date.now()), lastChange = useRef(Date.now());
  const signature = `${value.stage}\0${value.fraction}\0${value.detail}`;
  useEffect(() => { lastChange.current = Date.now(); }, [signature]);
  useEffect(() => { const timer = window.setInterval(() => setClock(Date.now()), 1000); return () => window.clearInterval(timer); }, []);
  const percent = Math.max(0, Math.min(100, Math.round(value.fraction * 100)));
  const stageLabels: Record<PipelineProgress["stage"], string> = {
    preprocessing: "Read filtering and sample assignment", umi: "UMI family grouping", consensus: "Consensus calling",
    contamination: "Contamination checks", postprocessing: "Alignment and downstream filtering", collapse: "Haplotype collapse",
    tree: "Phylogeny inference", complete: "Analysis complete",
  };
  const elapsed = Math.max(0, Math.floor((clock - started.current) / 1000)), quiet = Math.max(0, Math.floor((clock - lastChange.current) / 1000));
  const assignments = value.sampleAssignments ?? [], maximum = Math.max(1, ...assignments.map((row) => row.reads));
  const optional = (["contamination", "postprocessing", "collapse", "tree"] as const).includes(value.stage as OptionalStageName)
    ? value.stage as OptionalStageName : undefined;
  return <section className="run-progress" aria-live="polite">
    <div><span>{stageLabels[value.stage]}</span><span className="progress-actions"><strong>{percent}%</strong>{optional && value.fraction < 1 && onSkip && <button type="button" className="skip-stage" disabled={skipping} onClick={() => onSkip(optional)}>{skipping ? "Skipping…" : "Skip this step"}</button>}</span></div>
    <progress max="100" value={percent} />
    <p>{value.detail}</p>
    <div className="working-heartbeat"><i /><span>{quiet >= 3 ? `Still working · last pipeline update ${quiet.toLocaleString()} s ago` : "Working"}</span><em>Elapsed {Math.floor(elapsed / 60)}:{String(elapsed % 60).padStart(2, "0")}</em></div>
    {value.stage === "preprocessing" && assignments.length > 0 && <section className="demux-live" aria-label="Live reads assigned to each sample">
      <header><strong>Live sample assignments</strong><span>{assignments.reduce((sum, row) => sum + row.reads, 0).toLocaleString()} reads assigned</span></header>
      <div className="demux-bars">{assignments.map((row) => <div className="demux-row" key={row.sample}>
        <span title={row.sample}>{row.sample}</span><div><i style={{ width: `${row.reads / maximum * 100}%` }} /></div><strong>{row.reads.toLocaleString()}</strong>
      </div>)}</div>
    </section>}
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
  const [deferContamination, setDeferContamination] = useState(false);
  const [deferPostprocessing, setDeferPostprocessing] = useState(false);
  const [deferCollapse, setDeferCollapse] = useState(false);
  const [deferPhylogeny, setDeferPhylogeny] = useState(false);
  const [spoolStorage, setSpoolStorage] = useState<SpoolStorage>(() =>
    typeof (window as DirectoryPickerWindow).showDirectoryPicker === "function" ? "external-directory" : "automatic");
  const [scratchDirectory, setScratchDirectory] = useState<ExternalScratchDirectoryHandle>();
  const [progress, setProgress] = useState<PipelineProgress>();
  const [result, setResult] = useState<ResultBundle>();
  const [error, setError] = useState("");
  const [running, setRunning] = useState(false);
  const [cancelling, setCancelling] = useState(false);
  const [skippingStage, setSkippingStage] = useState<OptionalStageName>();
  const [navigationBlocked, setNavigationBlocked] = useState(false);
  const workerRef = useRef<Worker | undefined>(undefined);
  const protectWorkRef = useRef(false), historyGuardRef = useRef(false), confirmedLeaveRef = useRef(false);
  const maxWorkers = useMemo(() => Math.max(1, navigator.hardwareConcurrency || 1), []);
  const externalScratchSupported = typeof (window as DirectoryPickerWindow).showDirectoryPicker === "function";
  const noDownsampling = draftConfig.parameters.maxReadsPerSample === 0;
  const slots = useMemo(() => referenceSlots(draftConfig), [draftConfig]);
  const hasWorkToProtect = running || Boolean(result || fastq || configFile || configName !== "new configuration"
    || Object.keys(referenceAssignments).length || unassignedReferences.length);

  useEffect(() => { protectWorkRef.current = hasWorkToProtect; }, [hasWorkToProtect]);

  useEffect(() => {
    if (!hasWorkToProtect) return;
    const warnBeforeLeaving = (event: BeforeUnloadEvent) => {
      if (confirmedLeaveRef.current) return;
      event.preventDefault(); event.returnValue = "";
    };
    window.addEventListener("beforeunload", warnBeforeLeaving);
    return () => window.removeEventListener("beforeunload", warnBeforeLeaving);
  }, [hasWorkToProtect]);

  useEffect(() => {
    const interceptHistoryDeparture = () => {
      if (confirmedLeaveRef.current || !protectWorkRef.current) return;
      const previous = window.history.state && typeof window.history.state === "object" ? window.history.state : {};
      window.history.pushState({ ...previous, webporpidNavigationGuard: true }, "", window.location.href);
      historyGuardRef.current = true; setNavigationBlocked(true);
    };
    window.addEventListener("popstate", interceptHistoryDeparture);
    return () => window.removeEventListener("popstate", interceptHistoryDeparture);
  }, []);

  useEffect(() => {
    if (hasWorkToProtect && !historyGuardRef.current) {
      const previous = window.history.state && typeof window.history.state === "object" ? window.history.state : {};
      window.history.pushState({ ...previous, webporpidNavigationGuard: true }, "", window.location.href);
      historyGuardRef.current = true;
    } else if (!hasWorkToProtect && historyGuardRef.current) {
      historyGuardRef.current = false; window.history.back();
    }
  }, [hasWorkToProtect]);

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

  async function chooseScratchDirectory() {
    setError("");
    try {
      const pickerWindow = window as DirectoryPickerWindow;
      if (!pickerWindow.showDirectoryPicker) throw new Error("Direct external scratch directories are not supported by this browser. Use current Chrome/Edge or porpid-cli.");
      const directory = await pickerWindow.showDirectoryPicker({ id: "webporpid-scratch", mode: "readwrite" });
      const permission = directory.requestPermission ? await directory.requestPermission({ mode: "readwrite" }) : "granted";
      if (permission !== "granted") throw new Error("Write access to the selected scratch directory was not granted.");
      setScratchDirectory(directory); setSpoolStorage("external-directory");
    } catch (cause) {
      if (cause instanceof DOMException && cause.name === "AbortError") return;
      setError(cause instanceof Error ? cause.message : String(cause));
    }
  }

  async function run() {
    if (!fastq) { setError("Choose a FASTQ or FASTQ.GZ input file."); return; }
    if (!configText.trim()) { setError("Choose or paste a PORPID YAML configuration."); return; }
    setError(""); setResult(undefined); setRunning(true); setCancelling(false); setSkippingStage(undefined);
    setProgress({ stage: "preprocessing", fraction: 0, detail: "Checking the configuration and preparing local workers",
      sampleAssignments: draftConfig.samples.map((sample) => ({ sample: sample.name, reads: 0 })) });
    try {
      if (spoolStorage === "external-directory") {
        if (!scratchDirectory) throw new Error("Choose a writable scratch directory before starting the external-disk run.");
        const permission = scratchDirectory.queryPermission ? await scratchDirectory.queryPermission({ mode: "readwrite" }) : "granted";
        if (permission !== "granted") throw new Error("The selected scratch directory is no longer writable. Choose it again to restore permission.");
      }
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
        if (event.data.type === "progress") {
          const nextProgress = event.data.progress; setProgress(nextProgress);
          setSkippingStage((current) => current && (nextProgress.stage !== current || nextProgress.fraction >= 1) ? undefined : current);
        }
        if (event.data.type === "result") { setResult(event.data.result); setRunning(false); setCancelling(false); setSkippingStage(undefined); worker.terminate(); workerRef.current = undefined; }
        if (event.data.type === "error") { setError(event.data.message); setRunning(false); setCancelling(false); setSkippingStage(undefined); worker.terminate(); workerRef.current = undefined; }
      };
      worker.onerror = (event) => { setError(event.message || "The pipeline worker failed."); setRunning(false); setCancelling(false); worker.terminate(); workerRef.current = undefined; };
      worker.postMessage({ type: "run", file: fastq, config, workers, deferContamination, deferPostprocessing, deferCollapse, deferPhylogeny, inputMappings, spoolStorage,
        scratchDirectory: spoolStorage === "external-directory" ? scratchDirectory : undefined });
    } catch (cause) { setError(cause instanceof Error ? cause.message : String(cause)); setRunning(false); setCancelling(false); }
  }

  function cancel() {
    workerRef.current?.postMessage({ type: "cancel" });
    setCancelling(true); setError("");
    setProgress((current) => current ? { ...current, detail: "Cancelling safely and removing temporary read files…" } : current);
  }

  function skipOptionalStage(stage: OptionalStageName) {
    workerRef.current?.postMessage({ type: "skip-stage", stage }); setSkippingStage(stage); setError("");
    setProgress((current) => current ? { ...current, detail: `Skipping ${stage.replaceAll("-", " ")} safely; completed upstream results will be retained…` } : current);
  }

  function confirmPageDeparture() {
    confirmedLeaveRef.current = true; setNavigationBlocked(false); window.history.go(-2);
    window.setTimeout(() => { confirmedLeaveRef.current = false; }, 1200);
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
    <header className="site-header"><a className="brand" href={import.meta.env.BASE_URL} aria-label="webPORPID home"><span className="brand-mark">wp</span><span>webPORPID<small>Nanopore &amp; PacBio analysis</small></span></a><div className="site-header-right"><nav><a href="#run">Run</a><a href="./methods.html">Methods</a><a href="https://github.com/MurrellGroup/webPORPID">GitHub</a></nav><span className="app-version" title={`webPORPID ${APP_VERSION}`}>v{APP_VERSION}</span></div></header>
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
          <section className={`scratch-storage${noDownsampling ? " no-downsampling" : ""}`}>
            <header><div><span className="section-kicker">Temporary read storage</span><h3>Choose where raw-read partitions live</h3><MethodLink topic="streaming" /></div>{noDownsampling && <strong>External scratch recommended</strong>}</header>
            <div className="scratch-options">
              <label className={spoolStorage === "automatic" ? "selected" : ""}>
                <input type="radio" name="spool-storage" checked={spoolStorage === "automatic"} onChange={() => setSpoolStorage("automatic")} />
                <span><strong>Automatic browser storage</strong><small>Fast OPFS with adaptive compaction. Best when Maximum reads / sample is positive.</small></span>
              </label>
              <label className={spoolStorage === "external-directory" ? "selected recommended" : ""}>
                <input type="radio" name="spool-storage" checked={spoolStorage === "external-directory"} disabled={!externalScratchSupported}
                  onChange={() => setSpoolStorage("external-directory")} />
                <span><strong>External scratch directory{noDownsampling ? " · recommended" : ""}</strong><small>Streams partitions directly to a directory you choose, outside browser-origin quota. Temporary files are deleted after consensus; after a browser crash, remove any webporpid-scratch-* folder manually.</small></span>
              </label>
            </div>
            <div className="scratch-directory-row">
              <button type="button" className="secondary" disabled={!externalScratchSupported || running} onClick={() => void chooseScratchDirectory()}>{scratchDirectory ? "Change scratch directory" : "Choose scratch directory"}</button>
              <span>{scratchDirectory ? `${scratchDirectory.name} selected` : externalScratchSupported ? "No external directory selected" : "Unavailable in this browser; use current Chrome/Edge or porpid-cli"}</span>
            </div>
            {noDownsampling && <p className="scratch-recommendation"><strong>No downsampling is enabled.</strong> Every demultiplexed sequence and quality string must survive until the global UMI model and consensus pass. Use an external scratch disk with ample free space; 256 spool partitions is recommended for bounded per-worker memory on very large runs.</p>}
          </section>
          <div className="run-bar"><label><span>CPU workers</span><input type="number" min="1" max={maxWorkers} value={workers} onChange={(event) => setWorkers(Math.max(1, Math.min(maxWorkers, Number(event.target.value) || 1)))} /><small>Maximum detected: {maxWorkers}</small></label><div className="defer-options"><strong>Optional stages after consensus</strong><label><input type="checkbox" checked={deferContamination} onChange={(event) => setDeferContamination(event.target.checked)} /><span>Defer contamination checks</span></label><label><input type="checkbox" checked={deferPostprocessing} onChange={(event) => setDeferPostprocessing(event.target.checked)} /><span>Defer alignment + downstream filtering</span></label><label><input type="checkbox" checked={deferCollapse} onChange={(event) => setDeferCollapse(event.target.checked)} /><span>Defer haplotype collapse</span></label><label><input type="checkbox" checked={deferPhylogeny} onChange={(event) => setDeferPhylogeny(event.target.checked)} /><span>Defer phylogeny inference</span></label><small>Contamination can be deferred or skipped without blocking later work; downstream outputs then retain every sequence at that gate and are labelled unfiltered. Collapse requires downstream filtering, and the default phylogeny requires collapse.</small></div><div>{running ? <button className="danger" type="button" disabled={cancelling} onClick={cancel}>{cancelling ? "Cancelling…" : "Cancel analysis"}</button> : <button className="primary run-button" type="button" onClick={() => void run()}>Run webPORPID</button>}</div></div>
          {progress && running && <StageProgress value={progress} onSkip={skipOptionalStage} skipping={skippingStage === progress.stage} />}
          {error && <div className="error-box" role="alert">{error}</div>}
        </section>
        <section className="method-strip" id="about"><article><span>01</span><h3>Stream &amp; demultiplex</h3><p>Gzip chunks are decoded incrementally. Read-quality, primer, orientation, sample-ID and BPB logic follows the nanopore branch.</p><MethodLink topic="streaming" label="Detailed preprocessing methods" /></article><article><span>02</span><h3>Group &amp; call consensus</h3><p>Sparse two-error offspring likelihoods, LDA decisions, heteroduplex QC, seeded alignment and minimum-agreement counting.</p><MethodLink topic="consensus" label="Detailed consensus methods" /></article><article><span>03</span><h3>Filter &amp; explore</h3><p>Run-aware contamination, panel and functional filters, APOBEC model, aligned variants, phylogeny and component exports.</p><MethodLink topic="contamination" label="Detailed downstream methods" /></article></section>
      </main>
    </> : <ResultsExplorer bundle={result} onSaveResults={saveResults} onBundleChange={setResult} />}
    {result && <button className="new-analysis" type="button" onClick={() => { setResult(undefined); setProgress(undefined); }}>← New or load another analysis</button>}
    {slotConflict && <div className="modal-backdrop" role="presentation"><section className="slot-conflict" role="dialog" aria-modal="true" aria-labelledby="slot-conflict-title"><span className="section-kicker">Filename conflict</span><h3 id="slot-conflict-title">This file matches a different YAML slot</h3><p><strong>{slotConflict.file.name}</strong> was dropped onto <strong>{slotConflict.target.expectedName}</strong>, but its filename matches <strong>{slotConflict.matching.expectedName}</strong>.</p><div><button type="button" className="primary" onClick={() => resolveSlotConflict("swap")}>Swap into filename-matching slot</button><button type="button" onClick={() => resolveSlotConflict("matching")}>Use matching slot</button><button type="button" onClick={() => resolveSlotConflict("chosen")}>Keep chosen slot</button><button type="button" onClick={() => setSlotConflict(undefined)}>Cancel</button></div></section></div>}
    {navigationBlocked && <div className="modal-backdrop" role="presentation"><section className="slot-conflict navigation-warning" role="dialog" aria-modal="true" aria-labelledby="navigation-warning-title"><span className="section-kicker">Unsaved local work</span><h3 id="navigation-warning-title">Leave webPORPID?</h3><p>The Back action was paused because leaving now would discard the active analysis, selected input files, or loaded results. Save the results file first if you need to return to this work.</p><div><button type="button" className="primary" autoFocus onClick={() => setNavigationBlocked(false)}>Stay on this page</button><button type="button" className="danger" onClick={confirmPageDeparture}>Leave anyway</button></div></section></div>}
    <footer className="site-footer"><span>webPORPID · standalone local analysis</span><span>Deterministic, inspectable processing and session files</span></footer>
  </div>;
}
