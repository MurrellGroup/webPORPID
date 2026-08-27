import type { ConsensusRecord, UmiFamily } from "../types";

export function FamilySizeHistogram({ families }: { families: UmiFamily[] }) {
  const bins = Array(12).fill(0); families.forEach((family) => bins[Math.min(11, Math.floor(Math.log2(Math.max(1, family.familySize))))]++);
  const maximum = Math.max(1, ...bins), width = 520, height = 180, plotHeight = 130;
  return <svg className="report-chart" viewBox={`0 0 ${width} ${height}`} role="img" aria-label="UMI family-size histogram">
    <line x1="38" y1="10" x2="38" y2="140" /><line x1="38" y1="140" x2="510" y2="140" />
    {bins.map((count, index) => { const bar = count / maximum * plotHeight; return <g key={index}><rect x={43 + index * 38} y={140 - bar} width="29" height={bar} /><text x={57 + index * 38} y="158">{index < 11 ? `2^${index}` : "2^11+"}</text>{count > 0 && <text className="bar-count" x={57 + index * 38} y={Math.max(12, 136 - bar)}>{count}</text>}</g>; })}
    <text className="axis-label" x="275" y="177">family size (log₂ bins)</text>
  </svg>;
}

const COLORS: Record<string, string> = { A: "#4b9b62", C: "#337fb4", G: "#d99f19", T: "#d75b4e" };
export function AgreementPositionPlot({ consensuses }: { consensuses: ConsensusRecord[] }) {
  const sites = consensuses.flatMap((record) => record.lowAgreementSites), maximum = Math.max(1, ...sites.map((site) => site.position));
  return <svg className="report-chart" viewBox="0 0 520 180" role="img" aria-label="Minimum-agreement positions">
    <line x1="38" y1="10" x2="38" y2="140" /><line x1="38" y1="140" x2="510" y2="140" />
    {[0.25, 0.5, 0.75, 1].map((value) => <g key={value}><line className="grid" x1="38" y1={140 - value * 120} x2="510" y2={140 - value * 120} /><text x="32" y={144 - value * 120}>{value}</text></g>)}
    {sites.slice(0, 20_000).map((site, index) => <circle key={index} cx={38 + site.position / maximum * 472} cy={140 - site.agreement * 120}
      r={Math.min(6, 1.5 + site.modalRunLength * 0.45)} fill={COLORS[site.modalReadBase] ?? "#777"} fillOpacity=".56"><title>{site.modalReadBase} run {site.modalRunLength}; position {site.position}; agreement {site.agreement}</title></circle>)}
    <text className="axis-label" x="275" y="177">consensus position</text>
  </svg>;
}
