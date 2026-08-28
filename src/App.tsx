import { useMemo, useRef, useState } from "react";
import { blankConfig, parseConfigYaml, resolveReferenceFiles, serializeConfigYaml } from "./config";
import { decodeResult, encodeResult, safeDatasetName } from "./result-file";
import type { PipelineConfig, PipelineProgress, ResultBundle } from "./types";
import { ResultsExplorer } from "./components/results-explorer";
import { ConfigForm } from "./components/config-form";

function download(name: string, data: string | Uint8Array, mime: string) {
  const body: BlobPart = typeof data === "string" ? data : data.slice().buffer as ArrayBuffer;
  const url = URL.createObjectURL(new Blob([body], { type: mime }));
  const anchor = document.createElement("a");
  anchor.href = url; anchor.download = name; anchor.click();
  setTimeout(() => URL.revokeObjectURL(url), 0);
}

function localFiles(files: File[]) {
  const map = new Map<string, () => Promise<string>>();
  for (const file of files) {
    const relative = (file as File & { webkitRelativePath?: string }).webkitRelativePath;
    map.set(relative || file.name, () => file.text());
  }
  return map;
}

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
  const [configName, setConfigName] = useState("new configuration");
  const [configText, setConfigText] = useState(() => serializeConfigYaml(blankConfig()));
  const [configIssue, setConfigIssue] = useState("");
  const [references, setReferences] = useState<File[]>([]);
  const [workers, setWorkers] = useState(defaultWorkers);
  const [progress, setProgress] = useState<PipelineProgress>();
  const [result, setResult] = useState<ResultBundle>();
  const [error, setError] = useState("");
  const [running, setRunning] = useState(false);
  const [cancelling, setCancelling] = useState(false);
  const workerRef = useRef<Worker | undefined>(undefined);
  const maxWorkers = useMemo(() => Math.max(1, navigator.hardwareConcurrency || 1), []);

  async function chooseConfig(file?: File) {
    if (!file) return;
    const source = await file.text();
    setConfigName(file.name); setConfigText(source); setError("");
    try { setDraftConfig(parseConfigYaml(source)); setConfigIssue(""); }
    catch (cause) { setConfigIssue(cause instanceof Error ? cause.message : String(cause)); }
  }

  function changeConfig(config: PipelineConfig) {
    setDraftConfig(config); setConfigName("edited configuration");
    try { setConfigText(serializeConfigYaml(config)); setConfigIssue(""); }
    catch (cause) { setConfigIssue(cause instanceof Error ? cause.message : String(cause)); }
  }

  function changeRawYaml(source: string) {
    setConfigText(source); setConfigName("edited configuration");
    try { setDraftConfig(parseConfigYaml(source)); setConfigIssue(""); }
    catch (cause) { setConfigIssue(cause instanceof Error ? cause.message : String(cause)); }
  }

  function startFromScratch() {
    const config = blankConfig();
    setDraftConfig(config); setConfigText(serializeConfigYaml(config)); setConfigName("new configuration"); setConfigIssue("");
  }

  async function run() {
    if (!fastq) { setError("Choose a FASTQ or FASTQ.GZ input file."); return; }
    if (!configText.trim()) { setError("Choose or paste a PORPID YAML configuration."); return; }
    setError(""); setResult(undefined); setRunning(true); setCancelling(false);
    try {
      const config = await resolveReferenceFiles(parseConfigYaml(configText), localFiles(references));
      const worker = new Worker(new URL("./pipeline-worker.ts", import.meta.url), { type: "module" });
      workerRef.current = worker;
      worker.onmessage = (event: MessageEvent<{ type: "progress"; progress: PipelineProgress } | { type: "result"; result: ResultBundle } | { type: "error"; message: string }>) => {
        if (event.data.type === "progress") setProgress(event.data.progress);
        if (event.data.type === "result") { setResult(event.data.result); setRunning(false); setCancelling(false); worker.terminate(); workerRef.current = undefined; }
        if (event.data.type === "error") { setError(event.data.message); setRunning(false); setCancelling(false); worker.terminate(); workerRef.current = undefined; }
      };
      worker.onerror = (event) => { setError(event.message || "The pipeline worker failed."); setRunning(false); setCancelling(false); worker.terminate(); workerRef.current = undefined; };
      worker.postMessage({ type: "run", file: fastq, config, workers });
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
          <div className="input-grid">
            <label className="file-card"><span className="step">1</span><strong>Compressed reads</strong><small>FASTQ or FASTQ.GZ; decompressed incrementally.</small><input type="file" accept=".fastq,.fq,.gz,application/gzip" onChange={(event) => setFastq(event.target.files?.[0])} /><em>{fastq ? `${fastq.name} · ${(fastq.size / 1048576).toFixed(1)} MiB` : "Choose reads"}</em></label>
            <label className="file-card"><span className="step">2</span><strong>PORPID configuration</strong><small>Original nanopore-branch YAML is accepted.</small><input type="file" accept=".yaml,.yml,text/yaml" onChange={(event) => void chooseConfig(event.target.files?.[0])} /><em>{configName || "Choose YAML"}</em></label>
            <label className="file-card"><span className="step">3</span><strong>Reference FASTAs</strong><small>Panel, contamination, and optional functional references.</small><input type="file" accept=".fasta,.fa,.fas,.fna,text/plain" multiple onChange={(event) => setReferences([...event.target.files ?? []])} /><em>{references.length ? `${references.length} files selected` : "Choose references"}</em></label>
          </div>
          <ConfigForm config={draftConfig} onChange={changeConfig} onReset={startFromScratch} onDownload={saveConfig} />
          <details className="config-editor"><summary><span>Raw YAML</span><small>{configIssue ? "Needs attention" : "Synchronized with the form"}</small></summary><textarea value={configText} onChange={(event) => changeRawYaml(event.target.value)} spellCheck={false} placeholder="dataset:\n  samples:\n    ..." /></details>
          {configIssue && <div className="config-warning" role="status"><strong>Configuration is not ready</strong><span>{configIssue}</span><small>Continue editing in the form or raw YAML. The analysis will not start until the YAML is valid.</small></div>}
          <div className="run-bar"><label><span>CPU workers</span><input type="number" min="1" max={maxWorkers} value={workers} onChange={(event) => setWorkers(Math.max(1, Math.min(maxWorkers, Number(event.target.value) || 1)))} /><small>Maximum detected: {maxWorkers}</small></label><div>{running ? <button className="danger" type="button" disabled={cancelling} onClick={cancel}>{cancelling ? "Cancelling…" : "Cancel analysis"}</button> : <button className="primary run-button" type="button" onClick={() => void run()}>Run webPORPID</button>}</div></div>
          {progress && running && <StageProgress value={progress} />}
          {error && <div className="error-box" role="alert">{error}</div>}
        </section>
        <section className="method-strip" id="about"><article><span>01</span><h3>Stream &amp; demultiplex</h3><p>Gzip chunks are decoded incrementally. Read-quality, primer, orientation, sample-ID and BPB logic follows the nanopore branch.</p></article><article><span>02</span><h3>Group &amp; call consensus</h3><p>Sparse two-error offspring likelihoods, LDA decisions, heteroduplex QC, seeded alignment and minimum-agreement counting.</p></article><article><span>03</span><h3>Filter &amp; explore</h3><p>Run-aware contamination, panel and functional filters, APOBEC model, aligned variants, phylogeny and component exports.</p></article></section>
      </main>
    </> : <ResultsExplorer bundle={result} onSaveResults={saveResults} onBundleChange={setResult} />}
    {result && <button className="new-analysis" type="button" onClick={() => { setResult(undefined); setProgress(undefined); }}>← New or load another analysis</button>}
    <footer className="site-footer"><span>webPORPID · standalone local analysis</span><span>PORPID nanopore parity target · deterministic audit trail</span></footer>
  </div>;
}
