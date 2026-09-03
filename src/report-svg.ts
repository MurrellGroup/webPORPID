import type { AlignmentVariant } from "./alignment-utils.ts";
import type { PostprocRecord, ResultBundle, UmiFamily } from "./types.ts";
import { classicalMds } from "./mds.ts";
import { staticTreeHighlighterSvg } from "./static-tree-highlighter.ts";

export interface ReportSvg { extension: string; text: string }

const CATEGORY_COLORS: Record<string, string> = {
  "family-size-reject": "#6ba3ca", "maybe-artefact": "#f29b4b", likely_real: "#45a849",
  "UMI_len != 8": "#df6767", "LDA-rejects": "#a988c7", "minag-reject": "#9e766d", heteroduplex: "#df92c6",
};
const BASE_COLORS: Record<string, string> = { A: "#d9342b", C: "#2c7fb8", G: "#f28e2b", T: "#49a64d" };
const xml = (value: unknown) => String(value).replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll(">", "&gt;").replaceAll('"', "&quot;");

function document(width: number, height: number, title: string, body: string) {
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img"><title>${xml(title)}</title><style>text{font-family:Arial,sans-serif;font-size:10px;fill:#19302b}.axis{stroke:#233832;stroke-width:1}.grid{stroke:#dde4e0;stroke-width:1}.label{font-size:11px;font-weight:600}</style><rect width="100%" height="100%" fill="#fffdf8"/>${body}</svg>`;
}

function empty(title: string, detail: string) {
  return document(720, 210, title, `<text x="30" y="45" class="label">${xml(title)}</text><text x="30" y="75">${xml(detail)}</text>`);
}

function hashFraction(value: string) {
  let hash = 2166136261;
  for (let index = 0; index < value.length; index++) hash = Math.imul(hash ^ value.charCodeAt(index), 16777619);
  return (hash >>> 0) / 0xffffffff;
}

function quantile(values: number[], probability: number) {
  if (!values.length) return 0;
  const sorted = [...values].sort((a, b) => a - b), position = Math.max(0, Math.min(1, probability)) * (sorted.length - 1);
  const lower = Math.floor(position), fraction = position - lower;
  return sorted[lower] * (1 - fraction) + sorted[Math.min(lower + 1, sorted.length - 1)] * fraction;
}

function category(family: UmiFamily, artefactCutoff: number, agreementThreshold: number) {
  if (family.disposition !== "likely_real") return family.disposition;
  if (family.familySize < artefactCutoff) return "maybe-artefact";
  if ((family.minimumAgreement ?? 1) < agreementThreshold) return "minag-reject";
  return "likely_real";
}

function legend(categories: string[], x = 525, y = 24) {
  return categories.map((name, index) => `<g transform="translate(${x},${y + index * 18})"><circle r="4" fill="${CATEGORY_COLORS[name] ?? "#777"}"/><text x="10" y="3">${xml(name)}</text></g>`).join("");
}

function umiJitter(bundle: ResultBundle, sample: string): ReportSvg {
  const families = bundle.umiFamilies.filter((row) => row.sample === sample), rows = families.filter((row) => row.disposition !== "BPB-rejects");
  if (!rows.length) return { extension: "umi-family-decisions.svg", text: empty("UMI family decisions", "No UMI families are available for this sample.") };
  const configured = bundle.config.samples.find((row) => row.name === sample), summary = bundle.summaries.find((row) => row.sample === sample);
  const agreement = configured?.agreementOverride ?? bundle.config.parameters.agreementThreshold;
  const outlier = configured?.outlierQuantileOverride ?? bundle.config.parameters.outlierQuantile;
  const artefact = summary?.artefactCutoff ?? 0, likelySizes = families.filter((row) => row.disposition === "likely_real").map((row) => row.familySize);
  const q = Math.ceil(quantile(likelySizes, outlier));
  const maximum = Math.max(1, ...rows.map((row) => Math.log2(Math.max(1, row.familySize))), Math.log2(Math.max(1, q)));
  const minimumUmi = Math.min(6.5, ...rows.map((row) => row.umi.length - .5)), maximumUmi = Math.max(9.5, ...rows.map((row) => row.umi.length + .5));
  const x = (value: number) => 45 + value / maximum * 455, y = (value: number) => 155 - (value - minimumUmi) / Math.max(1, maximumUmi - minimumUmi) * 137;
  const categories = [...new Set(rows.map((row) => category(row, artefact, agreement)))];
  let body = `<line class="axis" x1="45" y1="18" x2="45" y2="155"/><line class="axis" x1="45" y1="155" x2="500" y2="155"/>`;
  for (const value of [7, 8, 9].filter((tick) => tick >= minimumUmi && tick <= maximumUmi)) body += `<line class="grid" x1="45" x2="500" y1="${y(value)}" y2="${y(value)}"/><text x="36" y="${y(value) + 4}">${value}</text>`;
  for (let value = 0; value <= Math.floor(maximum); value++) body += `<line class="grid" x1="${x(value)}" x2="${x(value)}" y1="18" y2="155"/><text x="${x(value)}" y="171" text-anchor="middle">${value}</text>`;
  body += `<line x1="${x(Math.log2(Math.max(1, artefact)))}" x2="${x(Math.log2(Math.max(1, artefact)))}" y1="18" y2="155" stroke="#e41a1c" stroke-width="2"/><line x1="${x(Math.log2(Math.max(1, q)))}" x2="${x(Math.log2(Math.max(1, q)))}" y1="18" y2="155" stroke="#ff9800" stroke-width="2"/>`;
  for (const row of rows.slice(0, 30_000)) {
    const decision = category(row, artefact, agreement), jitter = (hashFraction(`${row.sampleIndex}\0${row.umi}`) - .5) * .35;
    body += `<circle cx="${x(Math.log2(Math.max(1, row.familySize)))}" cy="${y(row.umi.length + jitter)}" r="3.2" fill="${CATEGORY_COLORS[decision] ?? "#777"}" fill-opacity=".7"><title>${xml(`${decision}; family size ${row.familySize}; UMI length ${row.umi.length}`)}</title></circle>`;
  }
  body += legend(categories) + `<g transform="translate(525,160)"><line x2="20" stroke="#ff9800" stroke-width="2"/><text x="27" y="4">quantile threshold</text><line x2="20" y1="18" y2="18" stroke="#e41a1c" stroke-width="2"/><text x="27" y="22">artefact threshold</text></g><text x="272" y="202" text-anchor="middle">log₂ UMI family size</text><text transform="translate(13,95) rotate(-90)" text-anchor="middle">UMI length</text>`;
  return { extension: "umi-family-decisions.svg", text: document(720, 210, `${sample} UMI family decisions`, body) };
}

function artefactJitter(bundle: ResultBundle, sample: string): ReportSvg {
  const families = bundle.umiFamilies.filter((row) => row.sample === sample), configured = bundle.config.samples.find((row) => row.name === sample);
  if (!families.length) return { extension: "artefact-cutoff.svg", text: empty("Artefact cutoff", "No UMI families are available for this sample.") };
  const agreement = configured?.agreementOverride ?? bundle.config.parameters.agreementThreshold;
  const outlier = configured?.outlierQuantileOverride ?? bundle.config.parameters.outlierQuantile;
  const fraction = configured?.artefactFractionOverride ?? bundle.config.parameters.artefactFraction;
  const artefact = bundle.summaries.find((row) => row.sample === sample)?.artefactCutoff ?? 0;
  const q = Math.ceil(quantile(families.filter((row) => row.disposition === "likely_real").map((row) => row.familySize), outlier));
  const allowed = new Set(["family-size-reject", "maybe-artefact", "likely_real", "minag-reject"]);
  const rows = families.filter((row) => row.familySize <= q && allowed.has(category(row, artefact, agreement))), maximum = Math.max(1, q), x = (value: number) => 45 + value / maximum * 455;
  const categories = ["family-size-reject", "maybe-artefact", "likely_real", "minag-reject"].filter((name) => rows.some((row) => category(row, artefact, agreement) === name));
  let body = `<line class="axis" x1="45" y1="20" x2="45" y2="150"/><line class="axis" x1="45" y1="150" x2="500" y2="150"/>`;
  for (let index = 0; index <= 5; index++) { const value = maximum * index / 5; body += `<line x1="${x(value)}" x2="${x(value)}" y1="145" y2="155"/><text x="${x(value)}" y="167" text-anchor="middle">${Math.round(value)}</text>`; }
  for (let index = 0; index < 10; index++) { const value = q * (index + .5) / 10; body += `<line x1="${x(value)}" x2="${x(value)}" y1="20" y2="150" stroke="#287bb2" stroke-width="1.4"/>`; }
  body += `<line x1="${x(artefact)}" x2="${x(artefact)}" y1="20" y2="150" stroke="#e41a1c" stroke-width="2"/><line x1="${x(q)}" x2="${x(q)}" y1="20" y2="150" stroke="#ff9800" stroke-width="2"/>`;
  for (const row of rows.slice(0, 30_000)) { const decision = category(row, artefact, agreement), jitter = (hashFraction(`artefact\0${row.sampleIndex}\0${row.umi}`) - .5) * 92; body += `<circle cx="${x(row.familySize)}" cy="${84 + jitter}" r="3.1" fill="${CATEGORY_COLORS[decision] ?? "#777"}" fill-opacity=".72"><title>${xml(`${decision}; family size ${row.familySize}`)}</title></circle>`; }
  body += legend(categories) + `<text x="272" y="177" text-anchor="middle">UMI family size</text><text x="270" y="199" text-anchor="middle">${fraction.toFixed(2)} artefact fraction (${artefact}) at quantile ${outlier.toFixed(2)} (${q})</text>`;
  return { extension: "artefact-cutoff.svg", text: document(720, 210, `${sample} artefact cutoff`, body) };
}

function lowAgreement(bundle: ResultBundle, sample: string): ReportSvg {
  const configured = bundle.config.samples.find((row) => row.name === sample), threshold = configured?.agreementOverride ?? bundle.config.parameters.agreementThreshold;
  const minimumFamily = bundle.summaries.find((row) => row.sample === sample)?.artefactCutoff ?? 0;
  const sites = bundle.consensuses.filter((row) => row.sample === sample && row.familySize >= minimumFamily).flatMap((record) => {
    const eligible = record.lowAgreementSites.filter((site) => site.agreement <= .9), maximumRun = Math.max(0, ...eligible.map((site) => site.modalRunLength));
    return eligible.filter((site) => site.modalRunLength === maximumRun).slice(0, 1);
  });
  if (!sites.length) return { extension: "low-agreement-positions.svg", text: empty("Low-agreement positions", "No qualifying low-agreement positions are available.") };
  const maximumPosition = Math.max(1, ...sites.map((site) => site.position)), minimumAgreement = Math.min(threshold, ...sites.map((site) => site.agreement), .9), lower = Math.max(0, minimumAgreement - .04), upper = .92;
  const x = (value: number) => 45 + value / maximumPosition * 455, y = (value: number) => 150 - (value - lower) / Math.max(.01, upper - lower) * 125;
  let body = `<line class="axis" x1="45" y1="20" x2="45" y2="150"/><line class="axis" x1="45" y1="150" x2="500" y2="150"/><line x1="45" x2="500" y1="${y(threshold)}" y2="${y(threshold)}" stroke="#7e168c" stroke-width="2"/>`;
  for (let index = 0; index <= 5; index++) { const value = maximumPosition * index / 5; body += `<line x1="${x(value)}" x2="${x(value)}" y1="145" y2="155"/><text x="${x(value)}" y="168" text-anchor="middle">${Math.round(value)}</text>`; }
  for (let index = 0; index < 5; index++) { const value = lower + (upper - lower) * index / 4; body += `<line class="grid" x1="42" x2="500" y1="${y(value)}" y2="${y(value)}"/><text x="37" y="${y(value) + 3}" text-anchor="end">${value.toFixed(2)}</text>`; }
  for (const site of sites.slice(0, 20_000)) body += `<circle cx="${x(site.position)}" cy="${y(site.agreement)}" r="${Math.min(8, 1.5 + site.modalRunLength * .55)}" fill="${BASE_COLORS[site.modalReadBase] ?? "#777"}" fill-opacity=".68"><title>${xml(`${site.modalReadBase}${site.modalRunLength}; position ${site.position}; agreement ${site.agreement}`)}</title></circle>`;
  body += `<text x="272" y="191" text-anchor="middle">sequence position</text><text transform="translate(13,96) rotate(-90)" text-anchor="middle">minimum agreement</text>`;
  return { extension: "low-agreement-positions.svg", text: document(720, 210, `${sample} low-agreement positions`, body) };
}

function jet(value: number) {
  const x = Math.max(0, Math.min(1, (value - .5) / .5)), red = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * x - 3))), green = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * x - 2))), blue = Math.max(0, Math.min(1, 1.5 - Math.abs(4 * x - 1)));
  return `rgb(${Math.round(red * 255)},${Math.round(green * 255)},${Math.round(blue * 255)})`;
}

function mds(bundle: ResultBundle, sample: string): ReportSvg {
  const eligible = bundle.records.filter((row) => row.sample === sample && row.alignedNt && row.artefactPass && row.agreementPass && row.contaminationPass && row.panelPass);
  const rows = [...eligible].sort((left, right) => right.familySize - left.familySize).slice(0, 600) as Array<PostprocRecord & { alignedNt: string }>;
  if (!rows.length) return { extension: "mds-apobec.svg", text: empty("Classical MDS", "No retained family sequences are available for MDS.") };
  const coordinates = classicalMds(rows.map((row) => row.alignedNt)), xs = coordinates.map((row) => row[0]), ys = coordinates.map((row) => row[1]);
  const minX = Math.min(...xs), maxX = Math.max(...xs), minY = Math.min(...ys), maxY = Math.max(...ys), span = Math.max(maxX - minX, maxY - minY, 1e-9);
  const displayMinX = (minX + maxX - span) / 2, displayMinY = (minY + maxY - span) / 2;
  const x = (value: number) => 55 + (value - displayMinX) / span * 440, y = (value: number) => 170 - (value - displayMinY) / span * 140;
  let body = `<line class="axis" x1="55" y1="25" x2="55" y2="170"/><line class="axis" x1="55" y1="170" x2="500" y2="170"/>`;
  for (let index = 0; index < 5; index++) { const fraction = index / 4, xv = displayMinX + fraction * span, yv = displayMinY + fraction * span; body += `<line class="grid" x1="${x(xv)}" x2="${x(xv)}" y1="25" y2="174"/><text x="${x(xv)}" y="185" text-anchor="middle">${xv.toExponential(1)}</text><line class="grid" x1="51" x2="500" y1="${y(yv)}" y2="${y(yv)}"/><text x="47" y="${y(yv) + 3}" text-anchor="end">${yv.toExponential(1)}</text>`; }
  rows.forEach((row, index) => { body += `<circle cx="${x(coordinates[index][0])}" cy="${y(coordinates[index][1])}" r="${Math.max(2, Math.min(11, Math.sqrt(row.familySize) * .65))}" fill="${jet(row.apobec?.posteriorGaInflated ?? .5)}" fill-opacity=".75" stroke="#24332f" stroke-width=".3"><title>${xml(`${row.id}; family size ${row.familySize}; APOBEC probability ${row.apobec?.posteriorGaInflated ?? .5}`)}</title></circle>`; });
  body += `<text x="275" y="205" text-anchor="middle">MDS 1</text><text transform="translate(16,105) rotate(-90)" text-anchor="middle">MDS 2</text>`;
  return { extension: "mds-apobec.svg", text: document(720, 230, `${sample} classical MDS`, body) };
}

/** Static, editable plot exports generated from the exact stored sample data. */
export function reportPlotSvgs(bundle: ResultBundle, sample: string, includeStaticTreeHighlighters = true): ReportSvg[] {
  const plots = [umiJitter(bundle, sample), artefactJitter(bundle, sample), lowAgreement(bundle, sample), mds(bundle, sample)];
  if (includeStaticTreeHighlighters) for (const variant of ["collapsed", "uncollapsed", "functional"] as AlignmentVariant[]) {
    const label = variant === "collapsed" ? "collapsed" : variant === "uncollapsed" ? "uncollapsed" : "functional";
    plots.push({ extension: `${label}-tree-highlighter.svg`, text: staticTreeHighlighterSvg(bundle, sample, variant) });
  }
  return plots;
}
