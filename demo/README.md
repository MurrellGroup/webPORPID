This demo is fully synthetic. The sequence names, panels, and reads are neutral toy data generated only to exercise primer matching, UMI grouping, indel-tolerant consensus, downstream alignment, export, and result loading.

No domain-specific or original-project demo material is bundled here.

Run the demo from the repository root:

```bash
node cli/porpid-cli.mjs run demo/synthetic_reads.fastq --config demo/synthetic_config.yaml --output demo.synthetic.webporpid --workers 2
node cli/porpid-cli.mjs inspect demo.synthetic.webporpid
node cli/porpid-cli.mjs export demo.synthetic.webporpid --component consensus-fasta --sample sample_1 --output consensus.fasta
node cli/porpid-cli.mjs export demo.synthetic.webporpid --component trimmed-aa-fasta --sample sample_1 --output trimmed-aa.fasta
```

The run produces one accepted `AACCGGTT` UMI family. The expected trimmed amino-acid sequence is `MPWAIGPYVYDGQLTTDNRQFVSEK*`.

The browser accepts the same YAML, FASTQ, and three FASTA files. For large inputs, `.fastq.gz` is supported and streamed in both the CLI and browser pipeline.
