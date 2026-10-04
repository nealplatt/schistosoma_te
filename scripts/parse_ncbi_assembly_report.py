"""Extract candidate metadata from local NCBI Datasets assembly reports.

Accepts a JSON object with a `reports` array, a single AssemblyDataReport object,
or JSONL with one AssemblyDataReport per nonblank line. Accepts camelCase report
fields and their snake_case CLI equivalents; conflicting aliases are rejected.
Unknown fields are ignored, unexpected containers rejected.
No network access, filtering, selection, or deduplication is performed.

Output is an intermediate TSV, NOT the curated assemblies.tsv schema. NA means
missing (absent, null, or empty string). Text is preserved verbatim using standard
CSV quoting for tabs/newlines; read with csv.DictReader(delimiter="\\t").
Assembly methods remain free text, and organism.infraspecificNames is used
without falling back to potentially conflicting BioSample attributes.
Relationship fields use only explicit NCBI metadata. Difference flags are
serialized as true/false (missing stays NA); difference reasons preserve the
pairedAssembly.differences text without combining other description fields.

Mapping reference (not fetched by this script):
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-assembly/
Input SHA-256, resolved local path, and 1-based record number retain traceability.
The original input must be retained for unextracted metadata and source provenance.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import sys


# Output field -> exact NCBI field path. Counts/lengths are reported in bp/counts.
field_paths = {
    "assembly_name": ("assemblyInfo", "assemblyName"),
    "assembly_level": ("assemblyInfo", "assemblyLevel"),
    "release_date": ("assemblyInfo", "releaseDate"),
    "reported_species_name": ("organism", "organismName"),
    "taxon_id": ("organism", "taxId"),
    "strain": ("organism", "infraspecificNames", "strain"),
    "isolate": ("organism", "infraspecificNames", "isolate"),
    "assembly_span_bp": ("assemblyStats", "totalSequenceLength"),
    "ungapped_span_bp": ("assemblyStats", "totalUngappedLength"),
    "contig_count": ("assemblyStats", "numberOfContigs"),
    "scaffold_count": ("assemblyStats", "numberOfScaffolds"),
    "contig_n50_bp": ("assemblyStats", "contigN50"),
    "scaffold_n50_bp": ("assemblyStats", "scaffoldN50"),
    "sequencing_technology": ("assemblyInfo", "sequencingTech"),
    "assembly_method": ("assemblyInfo", "assemblyMethod"),
    "representation_type": ("assemblyInfo", "assemblyType"),
    "assembly_status": ("assemblyInfo", "assemblyStatus"),
    "current_accession": ("currentAccession",),
    "paired_accession": ("pairedAccession",),
    "refseq_genbank_are_different": ("assemblyInfo", "pairedAssembly", "refseqGenbankAreDifferent"),
    "refseq_genbank_difference_reason": ("assemblyInfo", "pairedAssembly", "differences"),
}
numeric_fields = {"taxon_id", "assembly_span_bp", "ungapped_span_bp", "contig_count",
                  "scaffold_count", "contig_n50_bp", "scaffold_n50_bp"}
columns = ("assembly_accession_version", "assembly_accession", "assembly_version",
           *field_paths, "input_path", "input_sha256", "input_record_number")


def report_objects(value):
    """Recognize supported shapes without recursively guessing field locations."""
    if not isinstance(value, dict):
        raise ValueError("expected a JSON object")
    if "reports" in value:
        reports = value["reports"]
        if not isinstance(reports, list):
            raise ValueError("reports must be an array")
    else:
        reports = [value]
    for report in reports:
        if not isinstance(report, dict) or not any(
            field in report for field in ("accession", "assemblyInfo", "assemblyStats",
                                          "assembly_info", "assembly_stats", "organism")
        ):
            raise ValueError("expected an NCBI AssemblyDataReport object")
    return reports


def get_value(report, path):
    value = report
    for field in path:
        if value is None:
            return None
        if not isinstance(value, dict):
            raise ValueError(f"{'.'.join(path)}: expected an object before {field}")
        alias = re.sub(r"(?<!^)(?=[A-Z])", "_", field).lower()
        if alias != field and field in value and alias in value and value[field] != value[alias]:
            raise ValueError(f"conflicting aliases: {field} and {alias}")
        value = value.get(field) if field in value else value.get(alias)
    return value


def map_report(report):
    row = {}
    accession = report.get("accession")
    if accession is None or accession == "":
        row.update(assembly_accession_version="NA", assembly_accession="NA", assembly_version="NA")
    else:
        if not isinstance(accession, str):
            raise ValueError("accession must be a string")
        match = re.fullmatch(r"(GC[AF]_[0-9]+)(?:\.([0-9]+))?", accession)
        if not match:
            raise ValueError(f"unexpected assembly accession: {accession!r}")
        row.update(assembly_accession_version=accession,
                   assembly_accession=match[1], assembly_version=match[2] or "NA")
    for field, path in field_paths.items():
        value = get_value(report, path)
        if value is None or value == "":
            row[field] = "NA"
        elif field == "refseq_genbank_are_different":
            if not isinstance(value, bool):
                raise ValueError(f"{'.'.join(path)} must be a boolean")
            row[field] = str(value).lower()
        elif field in numeric_fields:
            if isinstance(value, bool) or not isinstance(value, (str, int)) or not re.fullmatch(
                r"[0-9]+", str(value)
            ):
                raise ValueError(f"{'.'.join(path)} must be a nonnegative integer or digit string")
            row[field] = str(value)
        elif isinstance(value, str):
            row[field] = value
        else:
            raise ValueError(f"{'.'.join(path)} must be a string")
    return row


def parse_report(input_path):
    """Return mapped records in input order, retaining duplicate records."""
    input_path = Path(input_path)
    raw = input_path.read_bytes()
    text = raw.decode("utf-8")
    if not text.strip():
        raise ValueError("input is empty; use {\"reports\": []} for an empty report")
    try:
        objects = [json.loads(text)]
    except json.JSONDecodeError:
        objects = []
        for line_number, line in enumerate(text.splitlines(), 1):
            if line.strip():
                try:
                    objects.append(json.loads(line))
                except json.JSONDecodeError as error:
                    raise ValueError(f"invalid JSON/JSONL at line {line_number}: {error.msg}") from error
    reports = []
    for value in objects:
        reports.extend(report_objects(value))
    rows = []
    for number, report in enumerate(reports, 1):
        try:
            row = map_report(report)
        except ValueError as error:
            raise ValueError(f"record {number}: {error}") from error
        row.update(input_path=str(input_path.resolve()),
                   input_sha256=hashlib.sha256(raw).hexdigest(), input_record_number=str(number))
        rows.append(row)
    return rows


def write_candidates(input_path, output_path):
    """Validate before writing; refuse to overwrite files or curated inventory."""
    output_path = Path(output_path)
    inventory_path = Path(__file__).resolve().parents[1] / "metadata" / "genomes"
    discovery_path = inventory_path / "discovery"
    if (output_path.resolve().is_relative_to(inventory_path)
            and not output_path.resolve().is_relative_to(discovery_path)):
        raise ValueError("candidate output under metadata/genomes must be in discovery")
    if output_path.resolve() == Path(input_path).resolve():
        raise ValueError("output must differ from input")
    rows = parse_report(input_path)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", type=Path, help="local NCBI assembly-report JSON or JSONL")
    parser.add_argument("output", type=Path,
                        help="new candidate TSV outside curated inventory (discovery allowed); parent directory must exist")
    args = parser.parse_args()
    try:
        count = write_candidates(args.input, args.output)
    except (OSError, UnicodeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    print(f"Wrote {count} candidate records to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
