export type MethodTopic =
  | "streaming" | "filtering" | "demultiplexing" | "subsampling" | "umi" | "heteroduplex"
  | "consensus" | "agreement" | "contamination" | "optional-stages" | "artefact"
  | "panel" | "functional" | "collapse" | "phylogeny" | "results";

const TARGETS: Record<MethodTopic, string> = {
  streaming: "methods-preprocessing.html#streaming-storage",
  filtering: "methods-preprocessing.html#read-filtering",
  demultiplexing: "methods-preprocessing.html#demultiplexing",
  subsampling: "methods-preprocessing.html#subsampling",
  umi: "methods-preprocessing.html#umi-offspring",
  heteroduplex: "methods-consensus.html#heteroduplex",
  consensus: "methods-consensus.html#family-consensus",
  agreement: "methods-consensus.html#agreement",
  contamination: "methods-downstream.html#contamination",
  "optional-stages": "methods-downstream.html#optional-stages",
  artefact: "methods-downstream.html#artefact-agreement",
  panel: "methods-downstream.html#panel-filter",
  functional: "methods-downstream.html#functional-filter",
  collapse: "methods-downstream.html#collapse",
  phylogeny: "methods-downstream.html#alignment-phylogeny",
  results: "methods-downstream.html#results-file",
};

export function methodHref(topic: MethodTopic) { return `./${TARGETS[topic]}`; }

export function MethodLink({ topic, label = "Method" }: { topic: MethodTopic; label?: string }) {
  return <a className="method-link" href={methodHref(topic)}>{label} <span aria-hidden="true">↗</span></a>;
}
