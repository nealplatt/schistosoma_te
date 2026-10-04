"""Fetch Schistosoma assembly metadata only, using Conda's ncbi_datasets environment.

Preserves stdout bytes unchanged and records command, timestamps, CLI version,
exit status, stderr, and SHA-256 in a JSON sidecar. No sequence/annotation
downloads, selection, deduplication, or quality filtering. NCBI's `current`
query may additionally return previous versions paired with current assemblies;
all returned records are retained. Requires installed CLI; never installs it.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys

from parse_ncbi_assembly_report import write_candidates


def fetch_metadata(output_directory):
    prefix = ["conda", "run", "--no-capture-output", "-n", "ncbi_datasets", "datasets"]
    version_result = subprocess.run(prefix + ["version"], capture_output=True, check=True)
    version_text = version_result.stdout.decode("utf-8").strip()
    match = re.fullmatch(r"datasets version: ([0-9]+\.[0-9]+\.[0-9]+)", version_text)
    if not match:
        raise ValueError(f"unexpected datasets version output: {version_text!r}")
    version = match[1]
    command = prefix + ["summary", "genome", "taxon", "Schistosoma",
                        "--assembly-source", "all", "--assembly-version", "current",
                        "--limit", "all", "--as-json-lines"]
    started = datetime.now(timezone.utc)
    stem = f"schistosoma_{started.strftime('%Y-%m-%dT%H%M%S%fZ')}_datasets-{version}"
    output_directory = Path(output_directory)
    raw_directory = output_directory / "raw"
    candidate_directory = output_directory / "candidates"
    raw_directory.mkdir(parents=True, exist_ok=True)
    candidate_directory.mkdir(parents=True, exist_ok=True)
    raw_path = raw_directory / f"{stem}.jsonl"
    provenance_path = raw_directory / f"{stem}.provenance.json"
    candidate_path = candidate_directory / f"{stem}.tsv"
    print(f"Executing: {shlex.join(command)}", flush=True)
    result = subprocess.run(command, capture_output=True)
    with raw_path.open("xb") as handle:
        handle.write(result.stdout)
    provenance = {
        "command": shlex.join(command), "command_argv": command,
        "query_taxon": "Schistosoma", "assembly_version_query": "current",
        "retrieval_timestamp": started.isoformat(),
        "completed_timestamp": datetime.now(timezone.utc).isoformat(),
        "datasets_version": version, "raw_response_path": str(raw_path.resolve()),
        "raw_response_checksum_algorithm": "sha256",
        "raw_response_checksum": hashlib.sha256(result.stdout).hexdigest(),
        "returncode": result.returncode,
        "stderr": result.stderr.decode("utf-8", errors="replace"),
        "parser_path": str(Path(__file__).with_name("parse_ncbi_assembly_report.py").resolve()),
        "parser_sha256": hashlib.sha256(
            Path(__file__).with_name("parse_ncbi_assembly_report.py").read_bytes()).hexdigest(),
        "candidate_path": str(candidate_path.resolve()),
    }
    with provenance_path.open("x", encoding="utf-8") as handle:
        json.dump(provenance, handle, indent=2)
        handle.write("\n")
    if result.returncode:
        raise ValueError(f"NCBI command failed ({result.returncode}); retained raw output and provenance at {raw_path}")
    count = write_candidates(raw_path, candidate_path)
    print(f"Raw response: {raw_path}\nProvenance: {provenance_path}\nCandidates: {candidate_path}\nCandidate rows: {count}")
    return raw_path, candidate_path, provenance_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path,
                        default=Path(__file__).resolve().parents[1] /
                        "metadata/genomes/discovery/ncbi")
    args = parser.parse_args()
    try:
        fetch_metadata(args.output_directory)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
