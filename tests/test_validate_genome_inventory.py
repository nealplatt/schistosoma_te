"""Synthetic inventory tests; no biological records or resource downloads."""

import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest


repository = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "validate_genome_inventory", repository / "scripts" / "validate_genome_inventory.py"
)
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)


class InventoryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        for name in validator.tables:
            self.write(name)

    def write(self, name, rows=(), columns=None):
        columns = validator.tables[name]["columns"] if columns is None else columns
        with (self.directory / f"{name}.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            for row in rows:
                writer.writerow([row.get(field, "NA") for field in columns])

    def assembly(self, assembly_id="test_assembly", **values):
        row = dict(assembly_id=assembly_id, source_namespace="synthetic",
                   source_record_url="https://example.invalid/assembly",
                   reported_species_name="synthetic taxon")
        row.update(values)
        return row

    def errors(self):
        return validator.validate_inventory(self.directory)

    def assert_error(self, text):
        self.assertTrue(any(text in error for error in self.errors()), self.errors())

    def test_valid_empty_inventory(self):
        self.assertEqual(self.errors(), [])

    def test_approved_headers(self):
        document = (repository / "docs" / "genome_inventory_schema.md").read_text()
        for name, table in validator.tables.items():
            section = document.split(f"## {name}.tsv\n", 1)[1].split("\n## ", 1)[0]
            fields = tuple(line.split("|")[1].strip().strip("`")
                           for line in section.splitlines() if line.startswith("| `"))
            self.assertEqual(table["columns"], fields)
            self.assertEqual((repository / "metadata" / "genomes" / f"{name}.tsv").read_text(),
                             "\t".join(fields) + "\n")

    def test_missing_required_column(self):
        self.write("assemblies", columns=validator.tables["assemblies"]["columns"][1:])
        self.assert_error("missing columns: assembly_id")

    def test_duplicate_assembly_id(self):
        self.write("assemblies", [self.assembly(), self.assembly()])
        self.assert_error("duplicate assembly_id")

    def test_broken_assembly_foreign_key(self):
        self.write("resources", [dict(resource_id="r", assembly_id="absent",
                                     resource_type="sequence", file_format="FASTA")])
        self.assert_error("assembly_id references unknown assembly")

    def test_negative_assembly_size(self):
        self.write("assemblies", [self.assembly(assembly_span_bp="-1")])
        self.assert_error("assembly_span_bp must be a nonnegative integer")

    def test_invalid_fraction(self):
        for field in validator.tables["assemblies"]["fractions"]:
            for value in ("-0.1", "1.1", "NaN", "Infinity", "not_numeric"):
                with self.subTest(field=field, value=value):
                    self.write("assemblies", [self.assembly(**{field: value})])
                    self.assert_error(f"{field} must be a finite fraction")

    def test_required_values(self):
        for value in ("", "NA", " "):
            with self.subTest(value=value):
                self.write("assemblies", [self.assembly(source_namespace=value)])
                self.assert_error("required field source_namespace is missing")

    def test_duplicate_external_identity(self):
        self.write("assemblies", [self.assembly("a", assembly_accession="X", assembly_version="1"),
                                  self.assembly("b", assembly_accession="X", assembly_version="1")])
        self.assert_error("duplicate external assembly identity")

    def test_missing_external_identities_are_not_duplicates(self):
        self.write("assemblies", [self.assembly("a"), self.assembly("b")])
        self.assertEqual(self.errors(), [])

    def test_broken_resource_references(self):
        self.write("assemblies", [self.assembly()])
        self.write("annotations", [dict(annotation_id="n", assembly_id="test_assembly",
                                        annotation_type="gene", resource_id="absent",
                                        compatibility_status="unchecked")])
        self.assert_error("resource_id references unknown resource")
        self.write("resources", [dict(resource_id="r", assembly_id="test_assembly",
                                     resource_type="sequence", file_format="FASTA",
                                     parent_resource_id="absent")])
        self.assert_error("parent_resource_id references unknown resource")

    def test_duplicate_selection_composite_key(self):
        self.write("assemblies", [self.assembly()])
        row = dict(analysis_id="analysis", comparison_unit="unit", assembly_id="test_assembly",
                   decision_status="undecided")
        self.write("assembly_selections", [row, row])
        self.assert_error("duplicate analysis_id/comparison_unit/assembly_id")

    def test_valid_populated_inventory(self):
        self.write("assemblies", [self.assembly(assembly_span_bp="0", gap_fraction="1",
                                               chromosome_assigned_fraction="0",
                                               quality_metric_scope="synthetic sequence",
                                               quality_metric_provenance="synthetic fixture")])
        self.write("resources", [dict(resource_id="r", assembly_id="test_assembly",
                                     resource_type="annotation", file_format="GFF")])
        self.write("annotations", [dict(annotation_id="n", assembly_id="test_assembly",
                                        annotation_type="gene", resource_id="r",
                                        compatibility_status="unchecked")])
        self.write("assembly_selections", [dict(analysis_id="analysis", comparison_unit="unit",
                                                assembly_id="test_assembly", decision_status="undecided")])
        self.assertEqual(self.errors(), [])

    def test_metric_provenance_is_conditionally_required(self):
        self.write("assemblies", [self.assembly(contig_count="1")])
        self.assert_error("conditionally required field quality_metric_provenance")

    def test_malformed_row(self):
        with (self.directory / "assemblies.tsv").open("a") as handle:
            handle.write("too\tfew\n")
        self.assert_error("row width does not match header")


if __name__ == "__main__":
    unittest.main()
