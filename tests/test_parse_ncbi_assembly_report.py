"""Offline tests using synthetic NCBI-shaped metadata only."""

import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import re


repository = Path(__file__).resolve().parents[1]
script = repository / "scripts" / "parse_ncbi_assembly_report.py"
fixtures = repository / "tests" / "fixtures" / "ncbi"
spec = importlib.util.spec_from_file_location("parse_ncbi_assembly_report", script)
parser = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parser)


class ParserTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    def test_exact_mapping_and_provenance(self):
        path = fixtures / "populated.json"
        rows = parser.parse_report(path)
        expected = {
            "assembly_accession_version": "GCA_999999991.2",
            "assembly_accession": "GCA_999999991", "assembly_version": "2",
            "assembly_name": "Synthetic assembly", "assembly_level": "Scaffold",
            "release_date": "2024-01-02", "reported_species_name": "Synthetic organism",
            "taxon_id": "999999999", "strain": "synthetic strain", "isolate": "synthetic isolate",
            "assembly_span_bp": "100000", "ungapped_span_bp": "99000",
            "contig_count": "20", "scaffold_count": "10",
            "contig_n50_bp": "10000", "scaffold_n50_bp": "20000",
            "sequencing_technology": "Synthetic short reads; synthetic long reads",
            "assembly_method": "Synthetic assembler v1; finishing method (unparsed)",
            "representation_type": "haploid-with-alt-loci",
            "assembly_status": "previous",
            "current_accession": "GCA_999999991.3",
            "paired_accession": "GCF_999999991.2",
            "refseq_genbank_are_different": "true",
            "refseq_genbank_difference_reason": "Synthetic difference: chromosome MT",
            "input_path": str(path.resolve()),
            "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "input_record_number": "1",
        }
        self.assertEqual(rows, [expected])
        self.assertEqual(set(parser.columns), set(expected))

    def test_substantial_missing_metadata(self):
        row = parser.parse_report(fixtures / "missing.json")[0]
        self.assertEqual(row["assembly_accession_version"], "GCF_999999992.1")
        self.assertEqual(row["assembly_accession"], "GCF_999999992")
        self.assertEqual(row["assembly_version"], "1")
        for field in parser.field_paths:
            self.assertEqual(row[field], "Synthetic missing metadata"
                             if field == "reported_species_name" else "NA")

    def test_snake_case_cli_fields(self):
        def snake(value):
            if isinstance(value, dict):
                return {re.sub(r"(?<!^)(?=[A-Z])", "_", k).lower(): snake(v)
                        for k, v in value.items()}
            if isinstance(value, list):
                return [snake(v) for v in value]
            return value
        original = json.loads((fixtures / "populated.json").read_text())["reports"][0]
        self.assertEqual(parser.map_report(snake(original)), parser.map_report(original))
        original["assembly_info"] = {"assembly_name": "conflicting"}
        with self.assertRaisesRegex(ValueError, "conflicting aliases"):
            parser.map_report(original)

    def test_explicit_difference_flag_false_and_missing(self):
        for value, expected in ((True, "true"), (False, "false"), (None, "NA")):
            with self.subTest(value=value):
                row = parser.map_report({"assembly_info": {"paired_assembly": {
                    "refseq_genbank_are_different": value}}})
                self.assertEqual(row["refseq_genbank_are_different"], expected)
                self.assertEqual(row["refseq_genbank_difference_reason"], "NA")
        for value in (0, 1, "true", [], {}):
            with self.subTest(value=value), self.assertRaises(ValueError):
                parser.map_report({"assembly_info": {"paired_assembly": {
                    "refseq_genbank_are_different": value}}})

    def test_relationships_are_not_inferred(self):
        row = parser.map_report({"accession": "GCA_999999991.2", "assemblyInfo": {
            "assemblyName": "Shared name", "pairedAssembly": {
                "accession": "GCF_999999991.2", "changed": "not a differences field"}}})
        for field in ("current_accession", "paired_accession", "assembly_status",
                      "refseq_genbank_are_different", "refseq_genbank_difference_reason"):
            self.assertEqual(row[field], "NA")

    def test_status_values_preserved_without_normalization(self):
        for value in ("current", "previous", "replaced", "suppressed", "other"):
            with self.subTest(value=value):
                self.assertEqual(parser.map_report({"assembly_info": {
                    "assembly_status": value}})["assembly_status"], value)

    def test_discovery_output_allowed_curated_output_blocked(self):
        fake_script = self.directory / "scripts" / "parser.py"
        output = self.directory / "metadata/genomes/discovery/ncbi/candidate.tsv"
        output.parent.mkdir(parents=True)
        with patch.object(parser, "__file__", str(fake_script)):
            self.assertEqual(parser.write_candidates(fixtures / "populated.json", output), 1)
            with self.assertRaises(ValueError):
                parser.write_candidates(fixtures / "populated.json",
                                        self.directory / "metadata/genomes/assemblies.tsv")

    def test_multiple_json_and_jsonl_preserve_duplicates_and_order(self):
        json_rows = parser.parse_report(fixtures / "multiple.json")
        jsonl_rows = parser.parse_report(fixtures / "multiple.jsonl")
        self.assertEqual(len(json_rows), 3)
        for rows in (json_rows, jsonl_rows):
            self.assertEqual([row["assembly_accession_version"] for row in rows],
                             ["GCA_999999991.2", "GCF_999999992.1", "GCA_999999991.2"])
            self.assertEqual([row["input_record_number"] for row in rows], ["1", "2", "3"])
        for a, b in zip(json_rows, jsonl_rows):
            self.assertEqual({k: v for k, v in a.items() if not k.startswith("input_")},
                             {k: v for k, v in b.items() if not k.startswith("input_")})

    def test_malformed_and_unexpected_fixtures(self):
        for filename in ("malformed.json", "unexpected.json"):
            with self.subTest(filename=filename), self.assertRaises(ValueError):
                parser.parse_report(fixtures / filename)

    def test_unexpected_shapes_and_types(self):
        path = self.directory / "input.json"
        for value in ({}, [], {"reports": {}}, {"reports": [None]},
                      {"organism": []}, {"accession": 42},
                      {"accession": "invalid"}, {"assemblyInfo": {"assemblyName": []}},
                      {"assemblyStats": {"numberOfContigs": True}},
                      {"assemblyStats": {"totalSequenceLength": -1}},
                      {"assemblyStats": {"contigN50": 1.5}}):
            with self.subTest(value=value):
                path.write_text(json.dumps(value))
                with self.assertRaises(ValueError):
                    parser.parse_report(path)

    def test_empty_report_and_missing_accession(self):
        path = self.directory / "input.json"
        path.write_text('{"reports": []}')
        self.assertEqual(parser.parse_report(path), [])
        path.write_text('{"organism": {"organismName": "Synthetic"}}')
        row = parser.parse_report(path)[0]
        for field in ("assembly_accession_version", "assembly_accession", "assembly_version"):
            self.assertEqual(row[field], "NA")
        path.write_text('{"accession": "GCA_999999993"}')
        self.assertEqual(parser.parse_report(path)[0]["assembly_version"], "NA")
        path.write_text("")
        with self.assertRaises(ValueError):
            parser.parse_report(path)

    def test_tsv_round_trip_free_text_and_zero(self):
        path = self.directory / "input.json"
        method = 'Synthetic method\twith "quotes"\nand a second line'
        path.write_text(json.dumps({"assemblyInfo": {"assemblyMethod": method},
                                    "assemblyStats": {"numberOfContigs": 0}}))
        output = self.directory / "candidates.tsv"
        self.assertEqual(parser.write_candidates(path, output), 1)
        with output.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(tuple(reader.fieldnames), parser.columns)
            row = next(reader)
            self.assertEqual(row["assembly_method"], method)
            self.assertEqual(row["contig_count"], "0")
            self.assertEqual(row["strain"], "NA")

    def test_no_biosample_or_current_accession_fallback(self):
        path = self.directory / "input.json"
        path.write_text(json.dumps({"organism": {}, "currentAccession": "GCA_999999993.2",
                                    "assemblyInfo": {"biosample": {"attributes": [
                                        {"name": "strain", "value": "do not infer"}]}}}))
        row = parser.parse_report(path)[0]
        self.assertEqual(row["strain"], "NA")
        self.assertEqual(row["assembly_accession_version"], "NA")

    def test_output_protection_and_invalid_input_leaves_no_output(self):
        output = self.directory / "candidates.tsv"
        with self.assertRaises(ValueError):
            parser.write_candidates(fixtures / "unexpected.json", output)
        self.assertFalse(output.exists())
        output.write_text("preserve me")
        with self.assertRaises(FileExistsError):
            parser.write_candidates(fixtures / "populated.json", output)
        self.assertEqual(output.read_text(), "preserve me")
        with self.assertRaises(ValueError):
            parser.write_candidates(fixtures / "populated.json",
                                    repository / "metadata" / "genomes" / "assemblies.tsv")

    def test_cli_help_success_and_failure(self):
        help_result = subprocess.run([sys.executable, "-B", str(script), "--help"],
                                     capture_output=True, text=True)
        self.assertEqual(help_result.returncode, 0)
        self.assertIn("intermediate TSV", help_result.stdout)
        output = self.directory / "candidates.tsv"
        success = subprocess.run([sys.executable, "-B", str(script),
                                  str(fixtures / "multiple.jsonl"), str(output)],
                                 capture_output=True, text=True)
        self.assertEqual(success.returncode, 0, success.stderr)
        self.assertIn("3 candidate records", success.stdout)
        for name in ("malformed.json", "unexpected.json", "absent.json"):
            invalid_output = self.directory / "invalid.tsv"
            result = subprocess.run([sys.executable, "-B", str(script),
                                     str(fixtures / name), str(invalid_output)],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("error:", result.stderr)
            self.assertFalse(invalid_output.exists())


if __name__ == "__main__":
    unittest.main()
