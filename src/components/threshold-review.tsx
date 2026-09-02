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

const DECISION_COLORS: Record<string, string> = {
  "family-size-reject": "#6ba3ca", "maybe-artefact": "#f29b4b", likely_real: "#45a849",
  "UMI_len != 8": "#df6767", "LDA-rejects": "#a988c7", "minag-reject": "#9e766d", heteroduplex: "#df92c6",
};
const DECISION_ORDER = ["family-size-reject", "maybe-artefact", "likely_real", "UMI_len != 8", "LDA-rejects", "minag-reject", "heteroduplex"];

function hashFraction(value: string) {
  let hash = 2166136261;
  for (let index = 0; index < value.length; index++) hash = Math.imul(hash ^ value.charCodeAt(index), 16777619);
  return (hash >>> 0) / 0xffffffff;
}

function liveUmiCategory(point: NonNullable<ThresholdReview["samples"][number]["displayPoints"]>[number], thresholds: Partial<Record<string, number>>) {
  if (point.disposition === "heteroduplex" || point.disposition === "BPB-rejects") return point.disposition;
  if ((point.posteriorProbability ?? 0) < (thresholds.ldaThreshold ?? 0)) return "LDA-rejects";
  if ((point.umi?.length ?? 8) !== 8) return "UMI_len != 8";
  if (point.familySize < (thresholds.familySizeThreshold ?? 0)) return "family-size-reject";
  return "likely_real";
}

function liveConsensusCategory(point: NonNullable<ThresholdReview["samples"][number]["displayPoints"]>[number], artefactCutoff: number, agreement: number) {
  if (point.disposition === "family-size-reject") return "family-size-reject";
  if (point.disposition !== "likely_real") return point.disposition;
  if (point.familySize < artefactCutoff) return "maybe-artefact";
  if ((point.minimumAgreement ?? 1) < agreement) return "minag-reject";
  return "likely_real";
}

function PlotLegend({ categories, x = 530, y = 25 }: { categories: string[]; x?: number; y?: number }) {
  return <g>{categories.map((category, index) => <g key={category} transform={`translate(${x},${y + index * 18})`}>
    <circle r="4" fill={DECISION_COLORS[category] ?? "#777"} /><text x="10" y="3">{category}</text>
  </g>)}</g>;
}

function DecisionPlot({ review, sample, thresholds }: { review: ThresholdReview; sample: ThresholdReview["samples"][number]; thresholds: Partial<Record<string, number>> }) {
  const points = sample.displayPoints ?? [], width = 720, height = 210, left = 45, plotRight = 500, top = 18, bottom = 155;
  if (review.phase === "umi") {
    const familyCutoff = thresholds.familySizeThreshold ?? 0;
    const maximumLog = Math.max(1, ...points.map((point) => Math.log2(Math.max(1, point.familySize))), Math.log2(Math.max(1, familyCutoff)));
    const minimumUmi = Math.min(6.5, ...points.map((point) => (point.umi?.length ?? 8) - .5));
    const maximumUmi = Math.max(9.5, ...points.map((point) => (point.umi?.length ?? 8) + .5));
    const x = (value: number) => left + value / maximumLog * (plotRight - left);
    const y = (value: number) => bottom - (value - minimumUmi) / Math.max(1, maximumUmi - minimumUmi) * (bottom - top);
    const categorized = points.map((point) => ({ point, category: liveUmiCategory(point, thresholds) }));
    const categories = DECISION_ORDER.filter((category) => categorized.some((row) => row.category === category));
    return <figure className="threshold-plot"><svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`Interactive UMI family-size decision plot for ${sample.sample}`}>
      <line x1={left} y1={top} x2={left} y2={bottom} /><line x1={left} y1={bottom} x2={plotRight} y2={bottom} />
      {[7, 8, 9].filter((value) => value >= minimumUmi && value <= maximumUmi).map((value) => <g key={value}><line x1={left} x2={plotRight} y1={y(value)} y2={y(value)} stroke="#e2e8e5" /><text x={left - 9} y={y(value) + 4} textAnchor="end">{value}</text></g>)}
      {Array.from({ length: Math.floor(maximumLog) + 1 }, (_, value) => <g key={value}><line x1={x(value)} x2={x(value)} y1={top} y2={bottom} stroke="#e2e8e5" /><text x={x(value)} y="171" textAnchor="middle">{value}</text></g>)}
      <line x1={x(Math.log2(Math.max(1, familyCutoff)))} x2={x(Math.log2(Math.max(1, familyCutoff)))} y1={top} y2={bottom} stroke="#e41a1c" strokeWidth="2" />
      {categorized.map(({ point, category }, index) => <circle key={`${point.umi ?? index}-${index}`} cx={x(Math.log2(Math.max(1, point.familySize)))}
        cy={y((point.umi?.length ?? 8) + (hashFraction(`${sample.sample}\0${point.umi ?? index}`) - .5) * .35)} r="3.2"
        fill={DECISION_COLORS[category] ?? "#777"} fillOpacity=".7"><title>{`${category}; family size ${point.familySize}; UMI length ${point.umi?.length ?? 8}; offspring probability ${(point.posteriorProbability ?? 0).toFixed(4)}`}</title></circle>)}
      <PlotLegend categories={categories} /><g transform="translate(530,160)"><line x1="0" x2="20" stroke="#e41a1c" strokeWidth="2" /><text x="27" y="4">family-size threshold {familyCutoff}</text></g>
      <text x="272" y="202" textAnchor="middle">log₂ UMI family size</text><text transform="translate(13,95) rotate(-90)" textAnchor="middle">UMI length</text>
      <text x="530" y="145">offspring probability ≥ {(thresholds.ldaThreshold ?? 0).toFixed(4)}</text>
    </svg><figcaption>{points.length.toLocaleString()} deterministic display points from {sample.totalFamilies.toLocaleString()} families. Colors update with the live thresholds and decisions are applied to the complete table. Heteroduplex status is computed during consensus, so it is not yet available at this pre-consensus checkpoint.</figcaption></figure>;
  }
  const quantileValue = Math.ceil(quantileFromCounts(sample.familySizeCounts, thresholds.outlierQuantile ?? 0));
  const artefactCutoff = Math.ceil(quantileFromCounts(sample.familySizeCounts, thresholds.outlierQuantile ?? 0) * (thresholds.artefactFraction ?? 0));
  const agreement = thresholds.agreementThreshold ?? 0, maximum = Math.max(1, quantileValue);
  const x = (value: number) => left + value / maximum * (plotRight - left);
  const allowed = new Set(["family-size-reject", "maybe-artefact", "likely_real", "minag-reject"]);
  const categorized = points.map((point) => ({ point, category: liveConsensusCategory(point, artefactCutoff, agreement) }))
    .filter(({ point, category }) => point.familySize <= quantileValue && allowed.has(category));
  const categories = DECISION_ORDER.filter((category) => categorized.some((row) => row.category === category));
  return <figure className="threshold-plot"><svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`Interactive artefact and minimum-agreement decision plot for ${sample.sample}`}>
    <line x1={left} y1="20" x2={left} y2="150" /><line x1={left} y1="150" x2={plotRight} y2="150" />
    {Array.from({ length: 6 }, (_, index) => maximum * index / 5).map((value, index) => <g key={index}><line x1={x(value)} x2={x(value)} y1="145" y2="155" /><text x={x(value)} y="167" textAnchor="middle">{Math.round(value)}</text></g>)}
    {Array.from({ length: 10 }, (_, index) => (index + .5) / 10).map((fraction) => <line key={fraction} x1={x(quantileValue * fraction)} x2={x(quantileValue * fraction)} y1="20" y2="150" stroke="#287bb2" strokeWidth="1.4" />)}
    <line x1={x(artefactCutoff)} x2={x(artefactCutoff)} y1="20" y2="150" stroke="#e41a1c" strokeWidth="2" /><line x1={x(quantileValue)} x2={x(quantileValue)} y1="20" y2="150" stroke="#ff9800" strokeWidth="2" />
    {categorized.map(({ point, category }, index) => <circle key={`${point.umi ?? index}-${index}`} cx={x(point.familySize)}
      cy={84 + (hashFraction(`artefact\0${sample.sample}\0${point.umi ?? index}`) - .5) * 92} r="3.1"
      fill={DECISION_COLORS[category] ?? "#777"} fillOpacity=".72"><title>{`${category}; family size ${point.familySize}; minimum agreement ${point.minimumAgreement ?? "not called"}`}</title></circle>)}
    <PlotLegend categories={categories} /><text x="272" y="177" textAnchor="middle">UMI family size</text>
    <text x="270" y="199" textAnchor="middle">artefact fraction {(thresholds.artefactFraction ?? 0).toFixed(3)} → cutoff {artefactCutoff}; quantile {(thresholds.outlierQuantile ?? 0).toFixed(3)} → {quantileValue}; min agreement {agreement.toFixed(3)}</text>
  </svg><figcaption>{categorized.length.toLocaleString()} deterministic display points. Family-size, artefact, and minimum-agreement classes update live; accepted thresholds are applied to every eligible consensus family.</figcaption></figure>;
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
