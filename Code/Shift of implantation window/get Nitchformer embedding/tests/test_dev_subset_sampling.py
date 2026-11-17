"""Development subset sampling tests."""

import numpy as np
import pandas as pd


def test_joint_sampling_is_reproducible_exact_and_ordered(load_script):
    module = load_script("01_create_dev_subset.py")
    obs = pd.DataFrame({
        "stage": np.repeat(["a", "b"], 50),
        "sample_id": np.tile(np.repeat(["s1", "s2"], 25), 2),
        "cell_type_major": np.tile(["stroma", "immune", "epithelial", "low_quality_cell"], 25),
    })
    first, manifest_a = module.select_indices(obs, 40, 20260806)
    second, manifest_b = module.select_indices(obs, 40, 20260806)
    np.testing.assert_array_equal(first, second)
    pd.testing.assert_frame_equal(manifest_a, manifest_b)
    assert len(first) == 40
    assert np.all(np.diff(first) > 0)
    selected = obs.iloc[first]
    assert set(selected["stage"]) == set(obs["stage"])
    assert set(selected["sample_id"]) == set(obs["sample_id"])
    assert set(selected["cell_type_major"]) == set(obs["cell_type_major"])


def test_rif_condition_is_derived_only_from_concordant_markers(load_script):
    module = load_script("01_create_dev_subset.py")
    obs = pd.DataFrame({
        "stage": ["pro_stage", "RIF_PO7"],
        "sample_id": ["pro_E1", "RIF_PO7_E2"],
    })
    derived = module.add_analysis_condition(obs)
    assert derived["analysis_condition"].tolist() == ["normal", "RIF"]

    bad = obs.copy()
    bad.loc[1, "sample_id"] = "sec_PO7_E2"
    with np.testing.assert_raises(ValueError):
        module.add_analysis_condition(bad)
