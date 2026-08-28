import { useMemo, useRef, type ReactNode } from "react";
import type { ConsensusRecord, PostprocRecord, UmiFamily } from "../types";

const CATEGORY_COLORS: Record<string, string> = {
  "family-size-reject": "#6ba3ca", "maybe-artefact": "#f29b4b", likely_real: "#45a849",
  "UMI_len != 8": "#df6767", "LDA-rejects": "#a988c7", "minag-reject": "#9e766d", heteroduplex: "#df92c6",
};
const BASE_COLORS: Record<string, string> = { A: "#d9342b", C: "#2c7fb8", G: "#f28e2b", T: "#49a64d" };

function ExportableSvg({ name, viewBox, label, children }: { name: string; viewBox: string; label: string; children: ReactNode }) {
  const ref = useRef<SVGSVGElement>(null);
  const save = () => {
    if (!ref.current) return;
    const clone = ref.current.cloneNode(true) as SVGSVGElement; clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(clone)], { type: "image/svg+xml;charset=utf-8" }));
    const anchor = document.createElement("a"); anchor.href = url; anchor.download = `${name}.svg`; anchor.click(); setTimeout(() => URL.revokeObjectURL(url), 0);
  };
  return <div className="plot-export"><button type="button" onClick={save}>Export SVG</button><svg ref={ref} className="report-chart" viewBox={viewBox} role="img" aria-label={label}>{children}</svg></div>;
}

function hashFraction(value: string): number {
  let hash = 2166136261;
  for (let index = 0; index < value.length; index += 1) hash = Math.imul(hash ^ value.charCodeAt(index), 16777619);
  return (hash >>> 0) / 0xffffffff;
}

function quantile(values: number[], probability: number): number {
  if (!values.length) return 0;
  const sorted = [...values].sort((a, b) => a - b), position = (sorted.length - 1) * probability;
  const lower = Math.floor(position), upper = Math.ceil(position), fraction = position - lower;
  return sorted[lower] * (1 - fraction) + sorted[upper] * fraction;
}

function familyCategory(family: UmiFamily, artefactCutoff: number, agreementThreshold: number): string {
  if (family.disposition !== "likely_real") return family.disposition;
  if (family.familySize < artefactCutoff) return "maybe-artefact";
  if ((family.minimumAgreement ?? 1) < agreementThreshold) return "minag-reject";
  return "likely_real";
}

function Legend({ categories, x = 525, y = 24 }: { categories: string[]; x?: number; y?: number }) {
  return <g>{categories.map((category, index) => <g key={category} transform={`translate(${x},${y + index * 18})`}><circle cx="0" cy="0" r="4" fill={CATEGORY_COLORS[category] ?? "#777"} /><text x="10" y="3">{category}</text></g>)}</g>;
}

export function UmiDecisionPlot({ families, artefactCutoff, agreementThreshold, outlierQuantile }: {
  families: UmiFamily[]; artefactCutoff: number; agreementThreshold: number; outlierQuantile: number;
}) {
  const rows = families.filter((family) => family.disposition !== "BPB-rejects"), likelySizes = families.filter((family) => family.disposition === "likely_real").map((family) => family.familySize);
  const quantileCutoff = Math.ceil(quantile(likelySizes, outlierQuantile)), maximum = Math.max(1, ...rows.map((row) => Math.log2(Math.max(1, row.familySize))), Math.log2(Math.max(1, quantileCutoff)));
  const categories = [...new Set(rows.map((row) => familyCategory(row, artefactCutoff, agreementThreshold)))];
  const minimumUmi = Math.min(6.5, ...rows.map((row) => row.umi.length - .5)), maximumUmi = Math.max(9.5, ...rows.map((row) => row.umi.length + .5));
  const x = (value: number) => 45 + value / Math.max(1, maximum) * 450, y = (value: number) => 155 - (value - minimumUmi) / Math.max(1, maximumUmi - minimumUmi) * 137;
  return <ExportableSvg name="umi-family-decisions" viewBox="0 0 720 210" label="UMI family size by UMI length and PORPID decision">
    <line x1="45" y1="18" x2="45" y2="155" /><line x1="45" y1="155" x2="500" y2="155" />
    {[7, 8, 9].filter((value) => value >= minimumUmi && value <= maximumUmi).map((value) => <g key={value}><line className="grid" x1="45" x2="500" y1={y(value)} y2={y(value)} /><text x="36" y={y(value) + 4}>{value}</text></g>)}
    {Array.from({ length: Math.floor(maximum) + 1 }, (_, value) => <g key={value}><line className="grid" x1={x(value)} x2={x(value)} y1="18" y2="155" /><text x={x(value)} y="171" textAnchor="middle">{value}</text></g>)}
    <line x1={x(Math.log2(Math.max(1, artefactCutoff)))} x2={x(Math.log2(Math.max(1, artefactCutoff)))} y1="18" y2="155" stroke="#e41a1c" strokeWidth="2" />
    <line x1={x(Math.log2(Math.max(1, quantileCutoff)))} x2={x(Math.log2(Math.max(1, quantileCutoff)))} y1="18" y2="155" stroke="#ff9800" strokeWidth="2" />
    {rows.slice(0, 30_000).map((row) => { const category = familyCategory(row, artefactCutoff, agreementThreshold), jitter = (hashFraction(`${row.sampleIndex}\0${row.umi}`) - .5) * .35; return <circle key={`${row.sampleIndex}-${row.umi}`} cx={x(Math.log2(Math.max(1, row.familySize)))} cy={y(row.umi.length + jitter)} r="3.2" fill={CATEGORY_COLORS[category] ?? "#777"} fillOpacity=".7"><title>{`${category}; family size ${row.familySize}; UMI length ${row.umi.length}`}</title></circle>; })}
    <Legend categories={categories} />
    <g transform="translate(525,160)"><line x1="0" x2="20" y1="0" y2="0" stroke="#ff9800" strokeWidth="2" /><text x="27" y="4">quantile threshold</text><line x1="0" x2="20" y1="18" y2="18" stroke="#e41a1c" strokeWidth="2" /><text x="27" y="22">artefact threshold</text></g>
    <text className="axis-label" x="272" y="202">log₂ UMI family size</text><text className="axis-label" transform="translate(13,95) rotate(-90)">UMI length</text>
  </ExportableSvg>;
}

export function ArtefactDecisionPlot({ families, artefactCutoff, agreementThreshold, outlierQuantile, artefactFraction }: {
  families: UmiFamily[]; artefactCutoff: number; agreementThreshold: number; outlierQuantile: number; artefactFraction: number;
}) {
  const likelySizes = families.filter((family) => family.disposition === "likely_real").map((family) => family.familySize), q = Math.ceil(quantile(likelySizes, outlierQuantile));
  const allowed = new Set(["family-size-reject", "maybe-artefact", "likely_real", "minag-reject"]);
  const rows = families.filter((family) => family.familySize <= q && allowed.has(familyCategory(family, artefactCutoff, agreementThreshold)));
  const maximum = Math.max(1, q), x = (value: number) => 45 + value / maximum * 455;
  const categories = ["family-size-reject", "maybe-artefact", "likely_real", "minag-reject"].filter((category) => rows.some((row) => familyCategory(row, artefactCutoff, agreementThreshold) === category));
  const maybe = rows.filter((row) => familyCategory(row, artefactCutoff, agreementThreshold) === "maybe-artefact").length;
  const real = rows.filter((row) => familyCategory(row, artefactCutoff, agreementThreshold) === "likely_real").length;
  const xTicks = Array.from({ length: 6 }, (_, index) => maximum * index / 5);
  return <ExportableSvg name="artefact-cutoff" viewBox="0 0 720 210" label="Family-size artefact cutoff jitter plot">
    <line x1="45" y1="20" x2="45" y2="150" /><line x1="45" y1="150" x2="500" y2="150" />
    {xTicks.map((value, index) => <g key={`tick-${index}`}><line className="grid" x1={x(value)} x2={x(value)} y1="145" y2="155" /><text x={x(value)} y="167" textAnchor="middle">{Math.round(value)}</text></g>)}
    {Array.from({ length: 10 }, (_, index) => (index + .5) / 10).map((fraction) => <line key={fraction} x1={x(q * fraction)} x2={x(q * fraction)} y1="20" y2="150" stroke="#287bb2" strokeWidth="1.4" />)}
    <line x1={x(artefactCutoff)} x2={x(artefactCutoff)} y1="20" y2="150" stroke="#e41a1c" strokeWidth="2" /><line x1={x(q)} x2={x(q)} y1="20" y2="150" stroke="#ff9800" strokeWidth="2" />
    {rows.slice(0, 30_000).map((row) => { const category = familyCategory(row, artefactCutoff, agreementThreshold), jitter = (hashFraction(`artifact\0${row.sampleIndex}\0${row.umi}`) - .5) * 92; return <circle key={`${row.sampleIndex}-${row.umi}`} cx={x(row.familySize)} cy={84 + jitter} r="3.1" fill={CATEGORY_COLORS[category] ?? "#777"} fillOpacity=".72"><title>{`${category}; family size ${row.familySize}`}</title></circle>; })}
    <Legend categories={categories} />
    <text className="axis-label" x="272" y="177">UMI family size</text><text className="axis-label" x="270" y="199">{artefactFraction.toFixed(2)} af-thresh ({artefactCutoff}) at quantile {outlierQuantile.toFixed(2)} ({q}) → {maybe + real ? (100 * maybe / (maybe + real)).toFixed(1) : "0.0"}% artefacts</text>
  </ExportableSvg>;
}

export function AgreementPositionPlot({ consensuses, threshold, minimumFamilySize = 0 }: { consensuses: ConsensusRecord[]; threshold: number; minimumFamilySize?: number }) {
  const sites = consensuses.filter((record) => record.familySize >= minimumFamilySize).flatMap((record) => {
    const eligible = record.lowAgreementSites.filter((site) => site.agreement <= .9);
    const maximumRun = Math.max(0, ...eligible.map((site) => site.modalRunLength));
    return eligible.filter((site) => site.modalRunLength === maximumRun).slice(0, 1);
  });
  const maximumPosition = Math.max(1, ...sites.map((site) => site.position)), minimumAgreement = Math.min(threshold, ...sites.map((site) => site.agreement), .9);
  const lowerAgreement = Math.max(0, minimumAgreement - .04), upperAgreement = .92;
  const x = (value: number) => 45 + value / maximumPosition * 455, y = (value: number) => 150 - (value - lowerAgreement) / Math.max(.01, upperAgreement - lowerAgreement) * 125;
  const xTicks = Array.from({ length: 6 }, (_, index) => maximumPosition * index / 5), yTicks = Array.from({ length: 5 }, (_, index) => lowerAgreement + (upperAgreement - lowerAgreement) * index / 4);
  return <ExportableSvg name="low-agreement-positions" viewBox="0 0 720 210" label="Minimum-agreement positions">
    <line x1="45" y1="20" x2="45" y2="150" /><line x1="45" y1="150" x2="500" y2="150" /><line x1="45" y1={y(threshold)} x2="500" y2={y(threshold)} stroke="#7e168c" strokeWidth="2" />
    {xTicks.map((value, index) => <g key={`x-${index}`}><line className="grid" x1={x(value)} x2={x(value)} y1="145" y2="155" /><text x={x(value)} y="168" textAnchor="middle">{Math.round(value)}</text></g>)}
    {yTicks.map((value, index) => <g key={`y-${index}`}><line className="grid" x1="42" x2="500" y1={y(value)} y2={y(value)} /><text x="37" y={y(value) + 3} textAnchor="end">{value.toFixed(2)}</text></g>)}
    {sites.slice(0, 20_000).map((site, index) => <circle key={index} cx={x(site.position)} cy={y(site.agreement)} r={Math.min(8, 1.5 + site.modalRunLength * .55)} fill={BASE_COLORS[site.modalReadBase] ?? "#777"} fillOpacity=".68"><title>{`${site.modalReadBase}${site.modalRunLength}; position ${site.position}; agreement ${site.agreement}`}</title></circle>)}
    <g transform="translate(525,24)">{["C", "G", "T", "A"].map((base, index) => <g key={base} transform={`translate(0,${index * 20})`}><circle r="5" fill={BASE_COLORS[base]} /><text x="12" y="4">{base} run length</text></g>)}<g transform="translate(0,92)"><line x1="-5" x2="20" stroke="#7e168c" strokeWidth="2" /><text x="27" y="4">minag threshold</text></g></g>
    <text className="axis-label" x="272" y="191">sequence position</text><text className="axis-label" transform="translate(13,96) rotate(-90)">minag</text>
  </ExportableSvg>;
}

function dinucleotideMatrix(families: UmiFamily[], weighted: boolean): number[][] {
  const bases = "ACGT", matrix = Array.from({ length: 4 }, () => Array(4).fill(0)), rows = families.filter((family) => family.disposition !== "BPB-rejects");
  let total = 0;
  for (const family of rows) for (let index = 0; index + 1 < family.umi.length; index += 1) {
    const left = bases.indexOf(family.umi[index]), right = bases.indexOf(family.umi[index + 1]); if (left < 0 || right < 0) continue;
    const weight = weighted ? family.familySize : 1; matrix[left][right] += weight; total += weight;
  }
  return matrix.map((row) => row.map((value) => total ? value / total : 0));
}

export function DinucleotideHeatmaps({ families }: { families: UmiFamily[] }) {
  const matrices = useMemo(() => [dinucleotideMatrix(families, false), dinucleotideMatrix(families, true)], [families]), bases = [..."ACGT"], displayRows = [3, 2, 1, 0];
  const color = (value: number) => { const amount = Math.max(0, Math.min(1, value / .2)); return `rgb(${Math.round(255 - 70 * amount)},${Math.round(247 - 187 * amount)},${Math.round(208 - 202 * amount)})`; };
  return <ExportableSvg name="umi-dinucleotide-frequencies" viewBox="0 0 720 265" label="Unweighted and read-weighted UMI dinucleotide frequencies">
    <defs><linearGradient id="dinucleotide-frequency" x1="0" y1="1" x2="0" y2="0"><stop offset="0" stopColor={color(0)} /><stop offset="1" stopColor={color(.2)} /></linearGradient></defs>
    {matrices.map((matrix, panel) => { const origin = 45 + panel * 350; return <g key={panel}><text x={origin + 104} y="22" textAnchor="middle" className="chart-title">{panel ? "all UMIs" : "UMI families"}</text>{displayRows.map((left, displayRow) => matrix[left].map((value, right) => <g key={`${left}-${right}`}><rect x={origin + right * 52} y={40 + displayRow * 46} width="52" height="46" fill={color(value)} /><text x={origin + right * 52 + 26} y={68 + displayRow * 46} textAnchor="middle">{value.toFixed(2)}</text></g>))}{bases.map((base, index) => <text key={`x-${base}`} x={origin + index * 52 + 26} y="242" textAnchor="middle">{base}</text>)}{displayRows.map((left, displayRow) => <text key={`y-${left}`} x={origin - 14} y={68 + displayRow * 46} textAnchor="middle">{bases[left]}</text>)}<rect x={origin + 220} y="40" width="12" height="184" fill="url(#dinucleotide-frequency)" /><text x={origin + 238} y="47">0.20</text><text x={origin + 238} y="226">0.00</text></g>; })}
    <text className="axis-label" x="350" y="260" textAnchor="middle">second UMI base (columns); first UMI base (rows)</text>
  </ExportableSvg>;
}

function pairwiseDistance(left: string, right: string): number {
  let differences = 0, compared = 0;
  for (let index = 0; index < Math.min(left.length, right.length); index += 1) if (left[index] !== "-" && right[index] !== "-") { compared += 1; if (left[index] !== right[index]) differences += 1; }
  return compared ? differences / compared : 0;
}

export function classicalMds(sequences: string[]): Array<[number, number]> {
  const count = sequences.length; if (count < 2) return sequences.map(() => [0, 0]);
  const squared = Array.from({ length: count }, () => new Float64Array(count)), means = new Float64Array(count); let total = 0;
  for (let left = 0; left < count; left += 1) for (let right = left + 1; right < count; right += 1) {
    const value = pairwiseDistance(sequences[left], sequences[right]) ** 2; squared[left][right] = squared[right][left] = value; means[left] += value; means[right] += value; total += value * 2;
  }
  for (let index = 0; index < count; index += 1) means[index] /= count; total /= count * count;
  const multiply = (vector: Float64Array) => {
    const output = new Float64Array(count);
    for (let row = 0; row < count; row += 1) for (let column = 0; column < count; column += 1)
      output[row] += -.5 * (squared[row][column] - means[row] - means[column] + total) * vector[column];
    return output;
  };
  let spectralShift = 0;
  for (let row = 0; row < count; row += 1) {
    let absoluteRow = 0;
    for (let column = 0; column < count; column += 1) absoluteRow += Math.abs(-.5 * (squared[row][column] - means[row] - means[column] + total));
    spectralShift = Math.max(spectralShift, absoluteRow);
  }
  const vectors: Float64Array[] = [], eigenvalues: number[] = [];
  for (let axis = 0; axis < 2; axis += 1) {
    let vector = Float64Array.from({ length: count }, (_, index) => Math.sin((index + 1) * (axis + 1) * 1.61803398875));
    for (let iteration = 0; iteration < 160; iteration += 1) {
      let mean = vector.reduce((sum, value) => sum + value, 0) / count; for (let index = 0; index < count; index += 1) vector[index] -= mean;
      let next = multiply(vector); for (let index = 0; index < count; index += 1) next[index] += spectralShift * vector[index];
      for (const prior of vectors) { let projection = 0; for (let index = 0; index < count; index += 1) projection += next[index] * prior[index]; for (let index = 0; index < count; index += 1) next[index] -= projection * prior[index]; }
      mean = next.reduce((sum, value) => sum + value, 0) / count; for (let index = 0; index < count; index += 1) next[index] -= mean;
      const norm = Math.hypot(...next); if (!norm) break; for (let index = 0; index < count; index += 1) next[index] /= norm;
      vector = next;
    }
    const product = multiply(vector); let eigenvalue = 0; for (let index = 0; index < count; index += 1) eigenvalue += vector[index] * product[index];
    vectors.push(vector); eigenvalues.push(Math.max(0, eigenvalue));
  }
  return sequences.map((_, index) => [vectors[0][index] * Math.sqrt(eigenvalues[0]), vectors[1][index] * Math.sqrt(eigenvalues[1])]);
}

function jet(value: number): string {
  const x = Math.max(0, Math.min(1, (value - .5) / .5)), red = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * x - 3))), green = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * x - 2))), blue = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * x - 1)));
  return `rgb(${Math.round(red * 255)},${Math.round(green * 255)},${Math.round(blue * 255)})`;
}

export function MdsApobecPlot({ records }: { records: PostprocRecord[] }) {
  const eligible = useMemo(() => records.filter((record) => record.alignedNt && record.artefactPass && record.agreementPass && record.contaminationPass && record.panelPass), [records]);
  const rows = useMemo(() => [...eligible].sort((left, right) => right.familySize - left.familySize).slice(0, 600), [eligible]);
  const coordinates = useMemo(() => classicalMds(rows.map((row) => row.alignedNt!)), [rows]);
  if (!rows.length) return <div className="empty-chart">No post-processing sequences passed for MDS.</div>;
  const xs = coordinates.map((row) => row[0]), ys = coordinates.map((row) => row[1]), minX = Math.min(...xs), maxX = Math.max(...xs), minY = Math.min(...ys), maxY = Math.max(...ys);
  const span = Math.max(maxX - minX, maxY - minY, 1e-9), displayMinX = (minX + maxX - span) / 2, displayMinY = (minY + maxY - span) / 2;
  const x = (value: number) => 55 + (value - displayMinX) / span * 440, y = (value: number) => 170 - (value - displayMinY) / span * 140;
  const ticks = Array.from({ length: 5 }, (_, index) => index / 4);
  return <ExportableSvg name="mds-apobec" viewBox="0 0 720 230" label="Classical multidimensional scaling colored by APOBEC probability">
    <defs><linearGradient id="jet-gradient" x1="0" y1="1" x2="0" y2="0"><stop offset="0" stopColor={jet(.5)} /><stop offset=".25" stopColor={jet(.625)} /><stop offset=".5" stopColor={jet(.75)} /><stop offset=".75" stopColor={jet(.875)} /><stop offset="1" stopColor={jet(1)} /></linearGradient></defs>
    <line x1="55" y1="25" x2="55" y2="170" /><line x1="55" y1="170" x2="500" y2="170" />
    {ticks.map((fraction) => { const value = displayMinX + fraction * span; return <g key={`x-${fraction}`}><line className="grid" x1={x(value)} x2={x(value)} y1="25" y2="174" /><text x={x(value)} y="185" textAnchor="middle">{value.toExponential(1)}</text></g>; })}
    {ticks.map((fraction) => { const value = displayMinY + fraction * span; return <g key={`y-${fraction}`}><line className="grid" x1="51" x2="500" y1={y(value)} y2={y(value)} /><text x="47" y={y(value) + 3} textAnchor="end">{value.toExponential(1)}</text></g>; })}
    {rows.map((row, index) => <circle key={row.id} cx={x(coordinates[index][0])} cy={y(coordinates[index][1])} r={Math.max(2, Math.min(11, Math.sqrt(row.familySize) * .65))} fill={jet(row.apobec?.posteriorGaInflated ?? .5)} fillOpacity=".75" stroke="#24332f" strokeWidth=".3"><title>{`family size ${row.familySize}; P(APOBEC|mutations) ${(row.apobec?.posteriorGaInflated ?? .5).toFixed(3)}`}</title></circle>)}
    <g transform="translate(532,25)"><text x="0" y="0">family size</text>{[10, 50, 100].map((value, index) => <g key={value} transform={`translate(12,${24 + index * 34})`}><circle r={Math.sqrt(value) * .65} fill="#555" fillOpacity=".7" /><text x="22" y="4">{value}</text></g>)}<rect x="104" y="12" width="18" height="125" fill="url(#jet-gradient)" /><text x="130" y="20">1.0</text><text x="130" y="78">0.75</text><text x="130" y="139">0.5</text></g>
    <text className="axis-label" x="275" y="205">MDS 1</text><text className="axis-label" transform="translate(16,105) rotate(-90)">MDS 2</text>{eligible.length > rows.length && <text x="55" y="222">Plot capped at the 600 largest families for interactive rendering.</text>}
  </ExportableSvg>;
}
