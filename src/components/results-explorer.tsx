import { useEffect, useMemo, useState } from "react";
import { exportComponent, type ExportKind, safeDatasetName } from "../result-file";
import type { ResultBundle } from "../types";
import { AgreementPositionPlot, FamilySizeHistogram } from "./charts";
import { AlignmentTreeViewer } from "./alignment-tree-viewer";

type Tab = "overview" | "families" | "sequences" | "contamination" | "alignment" | "log";
const PAGE_SIZE = 250;

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
    ["protein-alignment", "Protein alignment"], ["newick", "Newick tree"], ["log", "Run log"]];
  return <div className="export-menu"><select value={kind} onChange={(event) => setKind(event.target.value as ExportKind)}>{labels.map(([value, label]) => <option value={value} key={value}>{label}</option>)}</select>
    <button type="button" onClick={() => { const result = exportComponent(bundle, kind, sample); downloadText(`${safeDatasetName(bundle.config.dataset)}-${safeDatasetName(sample)}.${result.extension}`, result.text, result.mime); }}>Export</button></div>;
}

export function ResultsExplorer({ bundle, onSaveResults }: { bundle: ResultBundle; onSaveResults(): void }) {
  const [tab, setTab] = useState<Tab>("overview"), [sample, setSample] = useState(bundle.summaries[0]?.sample ?? bundle.config.samples[0]?.name ?? ""), [query, setQuery] = useState(""), [familyQuery, setFamilyQuery] = useState(""), [alphabet, setAlphabet] = useState<"nt" | "aa">("nt");
  const [familyPage, setFamilyPage] = useState(0), [sequencePage, setSequencePage] = useState(0), [contaminationPage, setContaminationPage] = useState(0);
  const summary = bundle.summaries.find((row) => row.sample === sample), families = useMemo(() => bundle.umiFamilies.filter((row) => row.sample === sample), [bundle, sample]), consensuses = useMemo(() => bundle.consensuses.filter((row) => row.sample === sample), [bundle, sample]);
  const filteredFamilies = useMemo(() => families.filter((row) => !familyQuery || row.umi.includes(familyQuery.toUpperCase()) || row.mostLikelyParent.includes(familyQuery.toUpperCase()) || row.disposition.toLowerCase().includes(familyQuery.toLowerCase())), [families, familyQuery]);
  const records = useMemo(() => bundle.records.filter((row) => row.sample === sample && (!query || row.id.toLowerCase().includes(query.toLowerCase()) || row.umi.includes(query.toUpperCase()))), [bundle, sample, query]);
  const contaminationRows = useMemo(() => bundle.contamination.filter((row) => row.sample === sample), [bundle, sample]);
  useEffect(() => { setFamilyPage(0); setSequencePage(0); setContaminationPage(0); }, [sample]);
  useEffect(() => setFamilyPage(0), [familyQuery]); useEffect(() => setSequencePage(0), [query]);
  const page = <T,>(rows: T[], index: number) => rows.slice(index * PAGE_SIZE, (index + 1) * PAGE_SIZE);
  const alignment = bundle.alignments[`${sample}/${alphabet === "nt" ? "nucleotide" : "protein"}`], tree = alphabet === "nt" ? bundle.trees[`${sample}/nucleotide`] : undefined;
  return <main className="results-page">
    <section className="results-hero"><div><span className="section-kicker">Loaded analysis</span><h1>{bundle.config.dataset}</h1><p>{bundle.provenance.inputName} · {bundle.provenance.createdUtc} · {bundle.provenance.workers} workers</p></div>
      <div className="result-actions"><select aria-label="Sample" value={sample} onChange={(event) => setSample(event.target.value)}>{bundle.summaries.map((row) => <option key={row.sample}>{row.sample}</option>)}</select><ExportMenu bundle={bundle} sample={sample} /><button className="primary" type="button" onClick={onSaveResults}>Save results file</button></div></section>
    <nav className="result-tabs">{(["overview", "families", "sequences", "contamination", "alignment", "log"] as Tab[]).map((value) => <button type="button" className={tab === value ? "active" : ""} onClick={() => setTab(value)} key={value}>{value}</button>)}</nav>
    {tab === "overview" && <section className="result-section">
      <div className="metric-grid"><article><span>Demultiplexed</span><strong>{summary?.demultiplexedReads.toLocaleString() ?? "0"}</strong></article><article><span>Observed UMIs</span><strong>{summary?.observedUmis.toLocaleString() ?? "0"}</strong></article><article><span>Consensus</span><strong>{summary?.consensusSequences.toLocaleString() ?? "0"}</strong></article><article><span>Postproc passed</span><strong>{summary?.postprocPassed.toLocaleString() ?? "0"}</strong></article><article><span>Artefact cutoff</span><strong>{summary?.artefactCutoff ?? 0}</strong></article></div>
      <Filters bundle={bundle} sample={sample} />
      <div className="chart-grid"><article><header><h3>Family-size distribution</h3><p>Log₂ bins retain the long tail without allowing extreme families to dominate.</p></header><FamilySizeHistogram families={families} /></article><article><header><h3>Low-agreement positions</h3><p>Point size is modal homopolymer run length; color is the modal read base.</p></header><AgreementPositionPlot consensuses={consensuses} /></article></div>
      <article className="provenance-card"><h3>Audit trail</h3><dl><div><dt>Input SHA-256</dt><dd>{bundle.provenance.inputSha256}</dd></div><div><dt>Config SHA-256</dt><dd>{bundle.provenance.configSha256}</dd></div><div><dt>Engine</dt><dd>{bundle.provenance.engine}</dd></div><div><dt>PORPID source</dt><dd>{bundle.provenance.upstreamBranch}@{bundle.provenance.upstreamCommit.slice(0, 12)}</dd></div></dl></article>
    </section>}
    {tab === "families" && <section className="result-section"><div className="table-heading"><div><h2>UMI family decisions</h2><p>Probabilistic offspring assignment, heteroduplex check, length and family-size gates.</p></div><input placeholder="Search UMI, parent, or class" value={familyQuery} onChange={(event) => setFamilyQuery(event.target.value)} /></div><div className="table-scroll"><table><thead><tr><th>UMI</th><th>fs</th><th>classification</th><th>parent</th><th>posterior</th><th>minag</th></tr></thead><tbody>{page(filteredFamilies, familyPage).map((row) => <tr key={`${row.sampleIndex}-${row.umi}`}><td><code>{row.umi}</code></td><td>{row.familySize}</td><td><span className={`status ${row.disposition === "likely_real" ? "pass" : "reject"}`}>{row.disposition}</span></td><td><code>{row.mostLikelyParent}</code></td><td>{row.posteriorProbability.toFixed(6)}</td><td>{row.minimumAgreement?.toFixed(2) ?? "—"}</td></tr>)}</tbody></table></div><Pager page={familyPage} count={filteredFamilies.length} onChange={setFamilyPage} /></section>}
    {tab === "sequences" && <section className="result-section"><div className="table-heading"><div><h2>Consensus and post-processing</h2><p>Every filter decision remains inspectable; FASTA exports use the exact stored sequences.</p></div><input placeholder="Search ID or UMI" value={query} onChange={(event) => setQuery(event.target.value)} /></div><div className="table-scroll"><table><thead><tr><th>Sequence</th><th>fs</th><th>minag</th><th>panel score</th><th>filters</th><th>APOBEC p</th><th>reason</th></tr></thead><tbody>{page(records, sequencePage).map((row) => <tr key={row.id}><td><code title={row.id}>{row.id}</code></td><td>{row.familySize}</td><td>{row.minimumAgreement.toFixed(2)}</td><td>{row.panelScore.toFixed(2)}</td><td><span className={`status ${row.artefactPass && row.agreementPass && row.contaminationPass && row.panelPass ? "pass" : "reject"}`}>{row.artefactPass && row.agreementPass && row.contaminationPass && row.panelPass ? "pass" : "reject"}</span></td><td>{row.apobec?.posteriorGaInflated.toFixed(3) ?? "—"}</td><td>{row.rejectionReasons.join("; ") || "—"}</td></tr>)}</tbody></table></div><Pager page={sequencePage} count={records.length} onChange={setSequencePage} /></section>}
    {tab === "contamination" && <section className="result-section"><div className="table-heading"><div><h2>Contamination report</h2><p>Primary calls and the wider zero-proportion suspect pass are retained separately.</p></div><span>{contaminationRows.length.toLocaleString()} calls</span></div><div className="table-scroll"><table><thead><tr><th>Sequence</th><th>Nearest non-self</th><th>distance</th><th>decision</th></tr></thead><tbody>{page(contaminationRows, contaminationPage).map((row, index) => <tr key={`${row.sequenceId}-${row.suspectOnly}-${contaminationPage}-${index}`}><td><code>{row.sequenceId}</code></td><td>{row.nearestNonselfVariant}</td><td>{row.nearestNonselfDistance.toPrecision(5)}</td><td><span className={`status ${row.discarded ? "reject" : "warn"}`}>{row.discarded ? "discarded" : row.suspectOnly ? "suspect pass" : "reported"}</span></td></tr>)}</tbody></table></div><Pager page={contaminationPage} count={contaminationRows.length} onChange={setContaminationPage} /></section>}
    {tab === "alignment" && <section className="result-section"><div className="alignment-switch"><div><h2>Tree + alignment workbench</h2><p>The tree and visible alignment rows share selection and vertical scrolling.</p></div><div><button className={alphabet === "nt" ? "active" : ""} onClick={() => setAlphabet("nt")}>Nucleotide</button><button className={alphabet === "aa" ? "active" : ""} onClick={() => setAlphabet("aa")}>Protein</button></div></div>{alignment ? <AlignmentTreeViewer fasta={alignment} newick={tree} alphabet={alphabet} /> : <div className="empty-state">No {alphabet === "nt" ? "nucleotide" : "protein"} alignment passed the selected sample’s filters.</div>}</section>}
    {tab === "log" && <section className="result-section"><div className="table-heading"><div><h2>Run log</h2><p>Persistent stage summaries and fallbacks stored inside the results file.</p></div></div><pre className="run-log">{bundle.log.join("\n")}</pre></section>}
  </main>;
}
