#!/usr/bin/env node
import { runCli } from "./porpid-cli.mjs";
runCli().catch((cause) => { process.stderr.write(`porpid-cli: ${cause instanceof Error ? cause.message : String(cause)}\n`); process.exitCode = 1; });
