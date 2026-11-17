"""Exact gene mapping tests."""

import pandas as pd


def test_model_lookup_strips_only_ensembl_versions(load_script):
    module = load_script("03_build_gene_mapping.py")
    var = pd.DataFrame({"gene_id": ["ENSG000001.4", ""], "gene_symbol": ["ABC1", "GENE-1"]})
    lookup = module.build_model_lookup(var, ["ENSG000001", "GENE-1"])
    assert lookup["ENSG000001"] == [0]
    assert lookup["GENE-1"] == [1]
    assert "gene" not in lookup


def test_ambiguous_model_identifiers_are_retained_as_ambiguous(load_script):
    module = load_script("03_build_gene_mapping.py")
    var = pd.DataFrame({"gene_symbol": ["DUP", "DUP"]})
    lookup = module.build_model_lookup(var, ["A", "B"])
    assert lookup["DUP"] == [0, 1]
