import type { PipelineConfig, PipelineParameters, SampleConfig } from "../types";
import { MethodLink, type MethodTopic } from "./method-link";

type NumberParameter = Exclude<keyof PipelineParameters, "contaminationFilter" | "deterministicSeed">;

interface NumericSetting {
  key: NumberParameter;
  label: string;
  help: string;
  min?: number;
  max?: number;
  step?: number;
}

const PREPROCESSING: NumericSetting[] = [
  { key: "errorRate", label: "FASTQ error rate", help: "Maximum expected read error rate.", min: 0, max: 0.99, step: 0.005 },
  { key: "minLength", label: "Minimum read length", help: "Reject shorter reads before demultiplexing.", min: 0, step: 1 },
  { key: "maxLength", label: "Maximum read length", help: "Reject longer reads before demultiplexing.", min: 1, step: 1 },
  { key: "primerTolerance", label: "Primer tolerance", help: "Allowed primer mismatches.", min: 0, step: 1 },
  { key: "primerWindow", label: "Primer search window", help: "Bases searched at each read end.", min: 16, step: 1 },
  { key: "primerChop", label: "Primer chop", help: "Additional bases removed beside primers.", min: 0, step: 1 },
  { key: "maxReadsPerSample", label: "Maximum reads / sample", help: "Deterministic cap after demultiplexing; zero keeps all reads and should use an external scratch directory for massive inputs.", min: 0, step: 1 },
];

const UMI_CONSENSUS: NumericSetting[] = [
  { key: "familySizeThreshold", label: "Family-size threshold", help: "Minimum reads in a UMI family.", min: 1, step: 1 },
  { key: "ldaThreshold", label: "Offspring probability", help: "LDA probability required to retain a UMI.", min: 0, max: 1, step: 0.001 },
  { key: "agreementThreshold", label: "Minimum agreement", help: "Minimum per-site consensus agreement.", min: 0, max: 1, step: 0.01 },
  { key: "artefactFraction", label: "Artefact fraction", help: "Within-sample abundance filter fraction.", min: 0, max: 1, step: 0.01 },
  { key: "outlierQuantile", label: "Outlier quantile", help: "Tail quantile used by the artefact model.", min: 0, max: 1, step: 0.001 },
];

const CONTAMINATION: NumericSetting[] = [
  { key: "contaminationClusterThreshold", label: "Cluster threshold", help: "Distance threshold for run clusters.", min: 0, max: 1, step: 0.001 },
  { key: "contaminationProportionThreshold", label: "Proportion threshold", help: "Non-self abundance threshold.", min: 0, max: 1, step: 0.01 },
  { key: "contaminationDistanceThreshold", label: "Distance threshold", help: "Nearest non-self sequence threshold.", min: 0, max: 1, step: 0.001 },
];

const POSTPROCESSING: NumericSetting[] = [
  { key: "panelThreshold", label: "Panel score threshold", help: "Maximum accepted misalignment score against the panel.", min: 0, step: 1 },
  { key: "functionalMatchThreshold", label: "Functional match", help: "Minimum translated reference match.", min: 0, max: 1, step: 0.01 },
];

function Field({ label, help, children, wide = false }: { label: string; help: string; children: React.ReactNode; wide?: boolean }) {
  return <label className={`config-field${wide ? " wide" : ""}`}><span>{label}</span>{children}<small>{help}</small></label>;
}

function NumericField({ setting, value, onChange }: { setting: NumericSetting; value: number; onChange: (value: number) => void }) {
  return <Field label={setting.label} help={setting.help}><input type="number" value={value} min={setting.min} max={setting.max} step={setting.step ?? "any"} onChange={(event) => {
    const next = Number(event.target.value); if (Number.isFinite(next)) onChange(next);
  }} /></Field>;
}

function OptionalNumber({ label, help, value, onChange, min, max, step = 0.01 }: {
  label: string; help: string; value?: number; onChange: (value: number | undefined) => void; min?: number; max?: number; step?: number;
}) {
  return <Field label={label} help={help}><input type="number" value={value ?? ""} placeholder="Use global" min={min} max={max} step={step} onChange={(event) => {
    const raw = event.target.value; onChange(raw === "" ? undefined : Number(raw));
  }} /></Field>;
}

function ParameterGroup({ title, description, settings, topic, config, onChange }: {
  title: string; description: string; settings: NumericSetting[]; topic: MethodTopic; config: PipelineConfig; onChange: (config: PipelineConfig) => void;
}) {
  const update = (key: NumberParameter, value: number) => onChange({ ...config, parameters: { ...config.parameters, [key]: value } });
  return <section className="parameter-group"><header><div><h4>{title}</h4><p>{description}</p></div><MethodLink topic={topic} /></header><div className="parameter-grid">
    {settings.map((setting) => <NumericField key={setting.key} setting={setting} value={config.parameters[setting.key]} onChange={(value) => update(setting.key, value)} />)}
  </div></section>;
}

function nextSampleName(samples: SampleConfig[]) {
  let index = samples.length + 1;
  while (samples.some((sample) => sample.name === `sample_${index}`)) index += 1;
  return `sample_${index}`;
}

export function ConfigForm({ config, onChange, onReset, onDownload }: {
  config: PipelineConfig; onChange: (config: PipelineConfig) => void; onReset: () => void; onDownload: () => void;
}) {
  const updateSample = <K extends keyof SampleConfig>(index: number, key: K, value: SampleConfig[K]) => {
    const samples = config.samples.map((sample, current) => current === index ? { ...sample, [key]: value } : sample);
    onChange({ ...config, samples });
  };
  const addSample = () => onChange({ ...config, samples: [...config.samples, {
    name: nextSampleName(config.samples), cdnaPrimer: "", secondStrandPrimer: "", panel: "panels/panel.fasta", panelSequences: [],
  }] });
  const removeSample = (index: number) => onChange({ ...config, samples: config.samples.filter((_, current) => current !== index) });

  return <div className="config-builder">
    <div className="config-builder-heading"><div><span className="section-kicker">Editable configuration</span><h3>Samples and pipeline settings</h3><p>Upload existing PORPID YAML or fill this in from scratch. The form and raw YAML stay synchronized.</p></div><div className="config-builder-actions"><button type="button" className="secondary" onClick={onDownload}>Download YAML</button><button type="button" className="secondary" onClick={onReset}>Start from scratch</button></div></div>
    <section className="parameter-group basics"><header><h4>Run identity &amp; references</h4><p>Paths are matched to the FASTA files selected above, including by filename.</p></header><div className="parameter-grid">
      <Field label="Dataset name" help="Used in the result filename and report." wide><input value={config.dataset} onChange={(event) => onChange({ ...config, dataset: event.target.value })} /></Field>
      <Field label="Contamination panel" help="FASTA used for external contamination checks." wide><input value={config.contaminationPanel} onChange={(event) => onChange({ ...config, contaminationPanel: event.target.value })} /></Field>
    </div></section>

    <section className="sample-section"><div className="sample-section-heading"><div><h4>Samples</h4><p>Lower-case bases in the cDNA primer identify the sample; <code>N</code> bases mark the UMI.</p></div><button type="button" onClick={addSample}>+ Add sample</button></div>
      <div className="sample-list">{config.samples.map((sample, index) => <article className="sample-card" key={index}>
        <header><span>{String(index + 1).padStart(2, "0")}</span><strong>{sample.name || "Unnamed sample"}</strong><button type="button" aria-label={`Remove ${sample.name || "sample"}`} disabled={config.samples.length === 1} onClick={() => removeSample(index)}>Remove</button></header>
        <div className="sample-fields">
          <Field label="Sample name" help="Unique label used throughout results."><input value={sample.name} onChange={(event) => updateSample(index, "name", event.target.value)} /></Field>
          <Field label="Panel FASTA" help="Reference alignment for this amplicon."><input value={sample.panel} onChange={(event) => updateSample(index, "panel", event.target.value)} /></Field>
          <Field label="cDNA primer" help="Preserve lower-case sample ID and N-marked UMI." wide><input className="sequence-input" value={sample.cdnaPrimer} onChange={(event) => updateSample(index, "cdnaPrimer", event.target.value)} spellCheck={false} /></Field>
          <Field label="Second-strand primer" help="Primer expected at the other read end." wide><input className="sequence-input" value={sample.secondStrandPrimer} onChange={(event) => updateSample(index, "secondStrandPrimer", event.target.value)} spellCheck={false} /></Field>
          <Field label="Functional reference" help="Optional coding reference FASTA."><input value={sample.functionalReference ?? ""} placeholder="Optional" onChange={(event) => updateSample(index, "functionalReference", event.target.value || undefined)} /></Field>
        </div>
        <details className="sample-overrides"><summary>Per-sample overrides</summary><div className="parameter-grid">
          <OptionalNumber label="Family size" help="Override the global family threshold." value={sample.familySizeOverride} min={1} step={1} onChange={(value) => updateSample(index, "familySizeOverride", value)} />
          <OptionalNumber label="Artefact fraction" help="Override the global artefact fraction." value={sample.artefactFractionOverride} min={0} max={1} onChange={(value) => updateSample(index, "artefactFractionOverride", value)} />
          <OptionalNumber label="Outlier quantile" help="Override the global outlier quantile." value={sample.outlierQuantileOverride} min={0} max={1} step={0.001} onChange={(value) => updateSample(index, "outlierQuantileOverride", value)} />
          <OptionalNumber label="Minimum agreement" help="Override consensus agreement." value={sample.agreementOverride} min={0} max={1} onChange={(value) => updateSample(index, "agreementOverride", value)} />
          <OptionalNumber label="Functional match" help="Override functional match threshold." value={sample.functionalMatchOverride} min={0} max={1} onChange={(value) => updateSample(index, "functionalMatchOverride", value)} />
        </div></details>
      </article>)}</div>
    </section>

    <div className="settings-groups">
      <ParameterGroup title="Preprocessing" description="Read quality, length, primer matching and downsampling." settings={PREPROCESSING} topic="filtering" config={config} onChange={onChange} />
      <ParameterGroup title="UMI &amp; consensus" description="Offspring classification, family consensus and artefact controls." settings={UMI_CONSENSUS} topic="umi" config={config} onChange={onChange} />
      <ParameterGroup title="Contamination" description="Run-aware self/non-self distance filters." settings={CONTAMINATION} topic="contamination" config={config} onChange={onChange} />
      <ParameterGroup title="Postprocessing" description="Reference-panel and optional coding-sequence filters." settings={POSTPROCESSING} topic="panel" config={config} onChange={onChange} />
      <section className="parameter-group"><header><h4>Resources &amp; reproducibility</h4><p>Disk partitioning and deterministic decisions.</p></header><div className="parameter-grid">
        <Field label="Spool partitions" help="Power of two; more partitions reduce peak memory."><select value={config.parameters.spoolPartitions} onChange={(event) => onChange({ ...config, parameters: { ...config.parameters, spoolPartitions: Number(event.target.value) } })}>{[1, 2, 4, 8, 16, 32, 64, 128, 256].map((value) => <option key={value}>{value}</option>)}</select></Field>
        <Field label="Deterministic seed" help="Unsigned integer stored in the audit trail."><input className="sequence-input" value={config.parameters.deterministicSeed.toString()} inputMode="numeric" onChange={(event) => { if (/^\d+$/.test(event.target.value)) onChange({ ...config, parameters: { ...config.parameters, deterministicSeed: BigInt(event.target.value) } }); }} /></Field>
        <label className="toggle-field"><input type="checkbox" checked={config.parameters.contaminationFilter} onChange={(event) => onChange({ ...config, parameters: { ...config.parameters, contaminationFilter: event.target.checked } })} /><span><strong>Run contamination filter</strong><small>When off, contamination classification is bypassed and every consensus is retained.</small></span></label>
      </div></section>
    </div>
  </div>;
}
