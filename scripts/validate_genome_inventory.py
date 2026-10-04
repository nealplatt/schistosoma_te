"""Validate inventory structure, not biological suitability.

Usage: python scripts/validate_genome_inventory.py [inventory_directory]
The default directory is metadata/genomes relative to the repository.
No resources are downloaded or modified. Full resource-content validation,
scientific review, and schema checks beyond this initial subset are deferred.
"""

import argparse
import csv
from decimal import Decimal, InvalidOperation
from pathlib import Path
import re


def is_missing(value):
    """NA is the documented sentinel; blank values are also absent."""
    return not value.strip() or value == "NA"


def validate_inventory(inventory_directory):
    """Return diagnostics; an empty list means the implemented checks pass."""
    inventory_directory = Path(inventory_directory)
    errors = []
    records = {}
    keys = {
        "assemblies": ("assembly_id",),
        "resources": ("resource_id",),
        "annotations": ("annotation_id",),
        "assembly_selections": ("analysis_id", "comparison_unit", "assembly_id"),
    }
    for name, spec in tables.items():
        records[name] = []
        path = inventory_directory / f"{name}.tsv"
        try:
            with path.open(encoding="utf-8", newline="") as handle:
                reader = csv.reader(handle, delimiter="\t")
                columns = next(reader, [])
                missing = set(spec["columns"]) - set(columns)
                if missing:
                    errors.append(f"{path.name}: missing columns: {', '.join(sorted(missing))}")
                if len(columns) != len(set(columns)):
                    errors.append(f"{path.name}: duplicate header columns")
                if missing or len(columns) != len(set(columns)):
                    continue
                seen = set()
                for values in reader:
                    location = f"{path.name}:{reader.line_num}"
                    if len(values) != len(columns):
                        errors.append(f"{location}: row width does not match header")
                        continue
                    row = dict(zip(columns, values))
                    records[name].append((location, row))
                    for field in spec["required"]:
                        if is_missing(row[field]):
                            errors.append(f"{location}: required field {field} is missing")
                    conditional = []
                    if name == "assemblies":
                        if any(not is_missing(row[field]) for field in
                               spec["integers"] + spec["fractions"]):
                            conditional.extend(("quality_metric_scope", "quality_metric_provenance"))
                        if not is_missing(row["related_assembly_id"]):
                            conditional.append("relationship_type")
                    elif name == "resources":
                        if any(not is_missing(row[field]) for field in
                               ("source_checksum", "local_checksum")):
                            conditional.append("checksum_algorithm")
                        if not is_missing(row["retrieval_date"]):
                            conditional.extend(("local_path", "local_checksum"))
                        if not is_missing(row["parent_resource_id"]):
                            conditional.extend(("transformation_provenance", "local_path", "local_checksum"))
                        if not is_missing(row["transformation_provenance"]):
                            conditional.append("parent_resource_id")
                    elif name == "annotations":
                        if not is_missing(row["target_assembly_accession"]):
                            conditional.append("target_assembly_namespace")
                    elif name == "assembly_selections":
                        if row["decision_status"] in ("proposed", "decided"):
                            conditional.append("rationale")
                        if row["decision_status"] == "decided":
                            conditional.append("decision_date")
                    for field in conditional:
                        if is_missing(row[field]):
                            errors.append(f"{location}: conditionally required field {field} is missing")
                    key = tuple(row[field] for field in keys[name])
                    if not any(is_missing(value) for value in key):
                        if key in seen:
                            errors.append(f"{location}: duplicate {'/'.join(keys[name])}: {key}")
                        seen.add(key)
                    for field in spec["integers"]:
                        value = row[field]
                        if not is_missing(value) and (
                            not re.fullmatch(r"[0-9]+", value)
                        ):
                            errors.append(f"{location}: {field} must be a nonnegative integer")
                    for field in spec["fractions"]:
                        value = row[field]
                        if is_missing(value):
                            continue
                        try:
                            number = Decimal(value)
                            valid = number.is_finite() and 0 <= number <= 1
                        except InvalidOperation:
                            valid = False
                        if not valid:
                            errors.append(f"{location}: {field} must be a finite fraction between 0 and 1")
        except (OSError, UnicodeError, csv.Error) as error:
            errors.append(f"{path.name}: cannot read TSV: {error}")

    assembly_ids = {row["assembly_id"] for _, row in records["assemblies"]
                    if not is_missing(row["assembly_id"])}
    resource_rows = {row["resource_id"]: row for _, row in records["resources"]
                     if not is_missing(row["resource_id"])}
    external_seen = set()
    for location, row in records["assemblies"]:
        identity = tuple(row[field] for field in
                         ("source_namespace", "assembly_accession", "assembly_version"))
        if not any(is_missing(value) for value in identity):
            if identity in external_seen:
                errors.append(f"{location}: duplicate external assembly identity: {identity}")
            external_seen.add(identity)

    for name, rows in records.items():
        for location, row in rows:
            assembly_fields = ("related_assembly_id",) if name == "assemblies" else ("assembly_id",)
            for field in assembly_fields:
                if not is_missing(row[field]) and row[field] not in assembly_ids:
                    errors.append(f"{location}: {field} references unknown assembly: {row[field]}")
            resource_fields = {"resources": ("parent_resource_id",),
                               "annotations": ("resource_id",)}.get(name, ())
            for field in resource_fields:
                value = row[field]
                if not is_missing(value) and value not in resource_rows:
                    errors.append(f"{location}: {field} references unknown resource: {value}")
            if name == "annotations" and row["resource_id"] in resource_rows:
                if resource_rows[row["resource_id"]]["assembly_id"] != row["assembly_id"]:
                    errors.append(f"{location}: annotation resource belongs to a different assembly")
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inventory_directory", nargs="?", type=Path,
                        default=Path(__file__).resolve().parents[1] / "metadata" / "genomes")
    args = parser.parse_args()
    errors = validate_inventory(args.inventory_directory)
    if errors:
        for error in errors:
            print(error)
        return 1
    print("Genome inventory validation passed.")
    return 0


# Approved columns in docs/genome_inventory_schema.md; kept explicit so the
# command works independently of Markdown parsing. Tests check schema agreement.
tables = {'assemblies': {'columns': ('assembly_id',
                            'source_namespace',
                            'assembly_accession',
                            'assembly_version',
                            'assembly_name',
                            'source_record_url',
                            'source_release',
                            'release_date',
                            'reported_species_name',
                            'accepted_species_name',
                            'taxon_id',
                            'taxonomy_source',
                            'strain',
                            'isolate',
                            'assembly_level',
                            'assembly_span_bp',
                            'ungapped_span_bp',
                            'contig_count',
                            'scaffold_count',
                            'contig_n50_bp',
                            'scaffold_n50_bp',
                            'chromosome_assigned_fraction',
                            'gap_fraction',
                            'quality_metric_scope',
                            'quality_metric_provenance',
                            'sequencing_technology',
                            'assembler',
                            'assembler_version',
                            'assembly_method_reference',
                            'representation_type',
                            'haplotype_label',
                            'related_assembly_id',
                            'relationship_type',
                            'notes'),
                'required': ('assembly_id',
                             'source_namespace',
                             'source_record_url',
                             'reported_species_name'),
                'integers': ('assembly_span_bp',
                             'ungapped_span_bp',
                             'contig_count',
                             'scaffold_count',
                             'contig_n50_bp',
                             'scaffold_n50_bp'),
                'fractions': ('chromosome_assigned_fraction', 'gap_fraction')},
 'resources': {'columns': ('resource_id',
                           'assembly_id',
                           'resource_type',
                           'sequence_component',
                           'masking_state',
                           'source_record_url',
                           'download_url',
                           'source_release',
                           'source_filename',
                           'local_path',
                           'retrieval_date',
                           'file_format',
                           'compression',
                           'file_size_bytes',
                           'checksum_algorithm',
                           'source_checksum',
                           'local_checksum',
                           'parent_resource_id',
                           'transformation_provenance',
                           'notes'),
               'required': ('resource_id',
                            'assembly_id',
                            'resource_type',
                            'file_format'),
               'integers': ('file_size_bytes',),
               'fractions': ()},
 'annotations': {'columns': ('annotation_id',
                             'assembly_id',
                             'annotation_type',
                             'provider',
                             'annotation_version',
                             'source_release',
                             'release_date',
                             'target_assembly_accession',
                             'target_assembly_version',
                             'target_assembly_namespace',
                             'resource_id',
                             'method_provenance',
                             'compatibility_status',
                             'notes'),
                 'required': ('annotation_id',
                              'assembly_id',
                              'annotation_type',
                              'resource_id',
                              'compatibility_status'),
                 'integers': (),
                 'fractions': ()},
 'assembly_selections': {'columns': ('analysis_id',
                                     'comparison_unit',
                                     'assembly_id',
                                     'role',
                                     'decision_status',
                                     'rationale',
                                     'decision_date'),
                         'required': ('analysis_id',
                                      'comparison_unit',
                                      'assembly_id',
                                      'decision_status'),
                         'integers': (),
                         'fractions': ()}}

if __name__ == "__main__":
    raise SystemExit(main())
