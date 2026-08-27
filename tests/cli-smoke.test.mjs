import assert from "node:assert/strict";
import { execFile } from "node:child_process";
import { mkdir, mkdtemp, readFile, rm, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";
import { promisify } from "node:util";
import test from "node:test";
import { gzipSync } from "node:zlib";

const execute = promisify(execFile);
const complements = { A: "T", C: "G", G: "C", T: "A", N: "N" };
const reverseComplement = (text) => [...text].reverse().map((base) => complements[base] ?? "N").join("");

test("porpid-cli runs, inspects, and exports a complete synthetic analysis", async () => {
  const testBase = resolve(".build/test-tmp"); await mkdir(testBase, { recursive: true });
  const directory = await mkdtemp(join(testBase, "webporpid-cli-test-"));
  try {
    await mkdir(join(directory, "panels"));
    const codons = ["GAC", "TTC", "AAG", "CCT", "GGA", "TAC", "CTG", "AAC", "GTC", "AGC", "TGG", "CCA", "GAT",
      "ACT", "CGT", "TCA", "AAA", "GCC", "TAT", "GGG", "ACC", "CAA", "GAG", "TCT", "ATC", "GTT", "CAG"];
    const payload = "ATG" + Array.from({ length: 40 }, (_, index) => codons[(index * 7 + 3) % codons.length]).join("") + "TAA";
    const payloads = [payload, payload, payload, `${payload.slice(0, 47)}A${payload.slice(47)}`, `${payload.slice(0, 81)}${payload.slice(82)}`];
    const fastq = payloads.map((variant, index) => {
      const oriented = `ACGTAAAACCGGTTGTCA${variant}`, sequence = `TAGG${reverseComplement(oriented)}${reverseComplement("CCGCT")}`;
      return `@synthetic_${index + 1}\n${sequence}\n+\n${"I".repeat(sequence.length)}\n`;
    }).join("");
    const reads = join(directory, "reads.fastq.gz"); await writeFile(reads, gzipSync(fastq));
    await writeFile(join(directory, "panels", "panel.fasta"), `>reference\n${reverseComplement(payload)}\n`);
    await writeFile(join(directory, "panels", "contam.fasta"), ">external\nACGTACGTACGTACGT\n");
    await writeFile(join(directory, "config.yaml"), `dataset: synthetic\nsamples:\n  sample_1:\n    cDNA_primer: CCGCTacgtaaNNNNNNNNGTCA\n    sec_str_primer: TAGG\n    panel: panels/panel.fasta\ncontaminationPanel: panels/contam.fasta\nparameters:\n  minLength: 20\n  maxLength: 300\n  primerTolerance: 0\n  primerWindow: 150\n  contaminationFilter: false\n  panelThreshold: 1000\n  spoolPartitions: 8\n`);
    const cli = resolve("cli/porpid-cli.mjs"), executable = process.env.WEBPORPID_CLI_EXECUTABLE || process.execPath;
    const prefix = process.env.WEBPORPID_CLI_EXECUTABLE ? [] : [cli], result = join(directory, "result.webporpid");
    await execute(executable, [...prefix, "run", reads, "--config", join(directory, "config.yaml"), "--output", result, "--workers", "2"], { timeout: 60_000 });
    const inspected = await execute(executable, [...prefix, "inspect", result]); const summary = JSON.parse(inspected.stdout);
    assert.equal(summary.quality.totalReads, 5); assert.equal(summary.quality.demultiplexedReads, 5);
    assert.equal(summary.components.consensuses, 1); assert.equal(summary.components.families, 1);
    const exported = join(directory, "consensus.fasta");
    await execute(executable, [...prefix, "export", result, "--component", "consensus-fasta", "--sample", "sample_1", "--output", exported]);
    const consensus = await readFile(exported, "utf8"); assert.match(consensus, /^>sample_1AACCGGTT fs=5 minag=/m);
    const serialResult = join(directory, "serial.webporpid"), serialExport = join(directory, "serial-consensus.fasta");
    await execute(executable, [...prefix, "run", reads, "--config", join(directory, "config.yaml"), "--output", serialResult, "--workers", "1"], { timeout: 60_000 });
    await execute(executable, [...prefix, "export", serialResult, "--component", "consensus-fasta", "--sample", "sample_1", "--output", serialExport]);
    assert.equal(await readFile(serialExport, "utf8"), consensus, "consensus changed with worker count");
  } finally { await rm(directory, { recursive: true, force: true }); }
});
