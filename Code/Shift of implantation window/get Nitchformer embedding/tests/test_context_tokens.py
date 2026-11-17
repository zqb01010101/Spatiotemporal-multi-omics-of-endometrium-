"""Fixed context policy tests."""

def test_no_assay_and_proxy_contexts(load_script):
    module = load_script("05_tokenize_nicheformer.py")
    assert module.CONTEXT["no_assay"] == ([5, 0, 4], [True, False, True])
    assert module.CONTEXT["10x3v3_proxy"] == ([5, 11, 4], [True, True, True])
    assert module.CONTEXT["no_assay"][0][0] == module.CONTEXT["10x3v3_proxy"][0][0]
    assert module.CONTEXT["no_assay"][0][2] == module.CONTEXT["10x3v3_proxy"][0][2]
