import { useMemo, useState } from "react";
import type { ThresholdReview, ThresholdSelection } from "../types";

type GlobalKey = "ldaThreshold" | "familySizeThreshold" | "artefactFraction" | "outlierQuantile" | "agreementThreshold";
type SampleKey = "familySizeOverride" | "artefactFractionOverride" | "outlierQuantileOverride" | "agreementOverride";
type Draft = { global: Partial<Record<GlobalKey, string>>; samples: Record<string, Partial<Record<SampleKey, string>>> };

function initialDraft(review: ThresholdReview): Draft {
  const global: Draft["global"] = {};
  for (const [key, value] of Object.entries(review.current) as Array<[GlobalKey, number | undefined]>) if (value !== undefined) global[key] = String(value);
  const samples: Draft["samples"] = {};
  for (const sample of review.samples) {
    samples[sample.sample] = review.phase === "umi"
      ? { familySizeOverride: sample.usesGlobal.familySizeThreshold ? "" : String(sample.current.familySizeThreshold) }
      : { artefactFractionOverride: sample.usesGlobal.artefactFraction ? "" : String(sample.current.artefactFraction),
        outlierQuantileOverride: sample.usesGlobal.outlierQuantile ? "" : String(sample.current.outlierQuantile),
        agreementOverride: sample.usesGlobal.agreementThreshold ? "" : String(sample.current.agreementThreshold) };
  }
  return { global, samples };
}

function parsed(value: string | undefined, fallback = 0) {
  const number = Number(value); return value !== undefined && value.trim() !== "" && Number.isFinite(number) ? number : fallback;
}

function NumericDecision({ label, help, raw, inherited, min, max, step, onChange }: {
  label: string; help: string; raw: string; inherited?: number; min: number; max: number; step: number; onChange(value: string): void;
}) {
  const effective = raw.trim() === "" ? inherited ?? min : parsed(raw, inherited ?? min);
  return <label className="threshold-control"><span><strong>{label}</strong><small>{help}</small></span><div>
    <input type="range" min={min} max={Math.max(min, max)} step={step} value={Math.max(min, Math.min(max, effective))}
      onChange={(event) => onChange(event.target.value)} />
    <input className="threshold-direct" type="text" inputMode="decimal" value={raw} placeholder={inherited === undefined ? "Value" : `Global (${inherited})`}
      aria-label={`${label} direct numeric value`} onChange={(event) => onChange(event.target.value)} />
  </div></label>;
}

function valueAt(counts: Array<[number, number]>, index: number) {
  let seen = 0;
  for (const [value, count] of counts) { if (index < seen + count) return value; seen += count; }
  return counts.at(-1)?.[0] ?? 0;
}

function quantileFromCounts(counts: Array<[number, number]>, probability: number) {
  const total = counts.reduce((sum, row) => sum + row[1], 0); if (!total) return 0;
  const position = Math.max(0, Math.min(1, probability)) * (total - 1), lower = Math.floor(position), fraction = position - lower;
  return valueAt(counts, lower) * (1 - fraction) + valueAt(counts, Math.min(lower + 1, total - 1)) * fraction;
}

function DecisionPlot({ review, sample, thresholds }: { review: ThresholdReview; sample: ThresholdReview["samples"][number]; thresholds: Partial<Record<string, number>> }) {
  const points = sample.displayPoints ?? [], width = 760, height = 300, left = 56, right = 18, top = 16, bottom = 44;
  const maximumFamily = Math.max(1, ...sample.familySizeCounts.map(([size]) => size));
  const x = (size: number) => left + Math.log1p(Math.max(0, Math.min(maximumFamily, size))) / Math.log1p(maximumFamily) * (width - left - right);
  const y = (value: number) => top + (1 - Math.max(0, Math.min(1, value))) * (height - top - bottom);
  const yValue = review.phase === "umi" ? "posteriorProbability" : "minimumAgreement";
  const horizontal = review.phase === "umi" ? thresholds.ldaThreshold ?? 0 : thresholds.agreementThreshold ?? 0;
  const vertical = review.phase === "umi" ? thresholds.familySizeThreshold ?? 0
    : Math.ceil(quantileFromCounts(sample.familySizeCounts, thresholds.outlierQuantile ?? 0) * (thresholds.artefactFraction ?? 0));
  return <figure className="threshold-plot"><svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`${review.phase === "umi" ? "UMI probability" : "Consensus agreement"} cutoff plot for ${sample.sample}`}>
    <rect x={left} y={top} width={width - left - right} height={height - top - bottom} fill="#fbfcfb" stroke="#c8d3cf" />
    {[0, .25, .5, .75, 1].map((tick) => <g key={tick}><line x1={left - 4} x2={width - right} y1={y(tick)} y2={y(tick)} stroke={tick ? "#e2e8e5" : "#879992"} /><text x={left - 9} y={y(tick) + 4} textAnchor="end">{tick}</text></g>)}
    {points.map((point, index) => <circle key={index} cx={x(point.familySize)} cy={y(Number(point[yValue]))} r="2.4" fill={point.disposition === "likely_real" ? "#08796f" : "#c5534f"} opacity=".45" />)}
    <line x1={left} x2={width - right} y1={y(horizontal)} y2={y(horizontal)} stroke="#b13c45" strokeWidth="2" strokeDasharray="7 5" />
    <line x1={x(vertical)} x2={x(vertical)} y1={top} y2={height - bottom} stroke="#4b3e91" strokeWidth="2" strokeDasharray="7 5" />
    <text x={(left + width - right) / 2} y={height - 8} textAnchor="middle">UMI family size (log scale)</text>
    <text transform={`translate(15 ${(top + height - bottom) / 2}) rotate(-90)`} textAnchor="middle">{review.phase === "umi" ? "Offspring probability" : "Minimum agreement"}</text>
    <text x={Math.min(width - 150, x(vertical) + 5)} y={top + 14} fill="#4b3e91">family cutoff {vertical}</text>
    <text x={left + 6} y={Math.max(top + 14, y(horizontal) - 6)} fill="#b13c45">{review.phase === "umi" ? "probability" : "agreement"} {horizontal}</text>
  </svg><figcaption>{points.length.toLocaleString()} display points from {sample.totalFamilies.toLocaleString()} families. Threshold decisions run over the complete family table.</figcaption></figure>;
}

export function ThresholdReviewDialog({ review, onAccept, onCancel }: {
  review: ThresholdReview; onAccept(selection: ThresholdSelection): void; onCancel(): void;
}) {
  const [draft, setDraft] = useState(() => initialDraft(review));
  const [sampleName, setSampleName] = useState(review.samples[0]?.sample ?? "");
  const [issue, setIssue] = useState("");
  const sample = review.samples.find((row) => row.sample === sampleName) ?? review.samples[0];
  const maximumFamily = Math.max(10, ...review.samples.flatMap((row) => row.familySizeCounts.map(([size]) => size)));
  const setGlobal = (key: GlobalKey, value: string) => setDraft((current) => ({ ...current, global: { ...current.global, [key]: value } }));
  const setSample = (key: SampleKey, value: string) => setDraft((current) => ({ ...current,
    samples: { ...current.samples, [sample.sample]: { ...current.samples[sample.sample], [key]: value } } }));
  const thresholds = useMemo(() => {
    const own = draft.samples[sample?.sample] ?? {};
    return review.phase === "umi" ? {
      ldaThreshold: parsed(draft.global.ldaThreshold),
      familySizeThreshold: parsed(own.familySizeOverride, parsed(draft.global.familySizeThreshold)),
    } : {
      artefactFraction: parsed(own.artefactFractionOverride, parsed(draft.global.artefactFraction)),
      outlierQuantile: parsed(own.outlierQuantileOverride, parsed(draft.global.outlierQuantile)),
      agreementThreshold: parsed(own.agreementOverride, parsed(draft.global.agreementThreshold)),
    };
  }, [draft, review.phase, sample?.sample]);

  const accept = () => {
    try {
      const number = (value: string | undefined, label: string, optional = false) => {
        if (optional && (value == null || !value.trim())) return undefined;
        const result = Number(value); if (!Number.isFinite(result)) throw new Error(`${label} must be a finite number.`); return result;
      };
      const parameters: ThresholdSelection["parameters"] = review.phase === "umi"
        ? { ldaThreshold: number(draft.global.ldaThreshold, "Offspring probability"),
          familySizeThreshold: number(draft.global.familySizeThreshold, "Global family size") }
        : { artefactFraction: number(draft.global.artefactFraction, "Artefact fraction"),
          outlierQuantile: number(draft.global.outlierQuantile, "Outlier quantile"),
          agreementThreshold: number(draft.global.agreementThreshold, "Minimum agreement") };
      if (parameters.familySizeThreshold !== undefined && (!Number.isSafeInteger(parameters.familySizeThreshold) || parameters.familySizeThreshold < 0))
        throw new Error("Family-size thresholds must be non-negative integers. Direct entry may exceed the slider range.");
      const samples = review.samples.map((row) => {
        const values = draft.samples[row.sample] ?? {};
        const output: ThresholdSelection["samples"][number] = { sample: row.sample };
        if (review.phase === "umi") output.familySizeOverride = number(values.familySizeOverride, `${row.sample} family size`, true);
        else {
          output.artefactFractionOverride = number(values.artefactFractionOverride, `${row.sample} artefact fraction`, true);
          output.outlierQuantileOverride = number(values.outlierQuantileOverride, `${row.sample} outlier quantile`, true);
          output.agreementOverride = number(values.agreementOverride, `${row.sample} minimum agreement`, true);
        }
        return output;
      });
      onAccept({ id: review.id, phase: review.phase, parameters, samples });
    } catch (cause) { setIssue(cause instanceof Error ? cause.message : String(cause)); }
  };

  if (!sample) return null;
  return <div className="modal-backdrop threshold-review-backdrop" role="presentation"><section className="threshold-review-dialog" role="dialog" aria-modal="true" aria-labelledby="threshold-review-title">
    <header><div><span className="section-kicker">Interactive decision checkpoint</span><h2 id="threshold-review-title">{review.title}</h2><p>{review.detail}</p></div><label><span>Plot sample</span><select value={sample.sample} onChange={(event) => setSampleName(event.target.value)}>{review.samples.map((row) => <option key={row.sample} value={row.sample}>{row.sample}{row.donorId ? ` · donor ${row.donorId}` : ""}</option>)}</select></label></header>
    <DecisionPlot review={review} sample={sample} thresholds={thresholds} />
    <div className="threshold-control-grid">
      {review.phase === "umi" ? <>
        <NumericDecision label="Global offspring probability" help="Posterior required before UMI-length and family-size checks." raw={draft.global.ldaThreshold ?? ""} min={0} max={1} step={.001} onChange={(value) => setGlobal("ldaThreshold", value)} />
        <NumericDecision label="Global family size" help="Minimum reads in a family; direct entry is not capped by this slider." raw={draft.global.familySizeThreshold ?? ""} min={0} max={maximumFamily} step={1} onChange={(value) => setGlobal("familySizeThreshold", value)} />
        <NumericDecision label={`${sample.sample} family size`} help="Leave direct entry blank to follow the global value." raw={draft.samples[sample.sample]?.familySizeOverride ?? ""} inherited={parsed(draft.global.familySizeThreshold)} min={0} max={maximumFamily} step={1} onChange={(value) => setSample("familySizeOverride", value)} />
      </> : <>
        <NumericDecision label="Global minimum agreement" help="Lowest accepted per-site family consensus agreement." raw={draft.global.agreementThreshold ?? ""} min={0} max={1} step={.01} onChange={(value) => setGlobal("agreementThreshold", value)} />
        <NumericDecision label="Global artefact fraction" help="Multiplier applied to the selected family-size quantile." raw={draft.global.artefactFraction ?? ""} min={0} max={1} step={.01} onChange={(value) => setGlobal("artefactFraction", value)} />
        <NumericDecision label="Global outlier quantile" help="Quantile of the contamination-eligible family-size distribution." raw={draft.global.outlierQuantile ?? ""} min={0} max={1} step={.001} onChange={(value) => setGlobal("outlierQuantile", value)} />
        <NumericDecision label={`${sample.sample} minimum agreement`} help="Leave blank to follow the global value." raw={draft.samples[sample.sample]?.agreementOverride ?? ""} inherited={parsed(draft.global.agreementThreshold)} min={0} max={1} step={.01} onChange={(value) => setSample("agreementOverride", value)} />
        <NumericDecision label={`${sample.sample} artefact fraction`} help="Leave blank to follow the global value." raw={draft.samples[sample.sample]?.artefactFractionOverride ?? ""} inherited={parsed(draft.global.artefactFraction)} min={0} max={1} step={.01} onChange={(value) => setSample("artefactFractionOverride", value)} />
        <NumericDecision label={`${sample.sample} outlier quantile`} help="Leave blank to follow the global value." raw={draft.samples[sample.sample]?.outlierQuantileOverride ?? ""} inherited={parsed(draft.global.outlierQuantile)} min={0} max={1} step={.001} onChange={(value) => setSample("outlierQuantileOverride", value)} />
      </>}
    </div>
    <p className="threshold-unbounded-note">Slider ranges are convenient guides only. Direct numeric fields have no slider-imposed limits; scientifically invalid non-numeric values are rejected.</p>
    {issue && <div className="error-box" role="alert">{issue}</div>}
    <footer><button type="button" className="danger" onClick={onCancel}>Cancel analysis</button><button type="button" className="primary" onClick={accept}>Accept thresholds and continue</button></footer>
  </section></div>;
}
