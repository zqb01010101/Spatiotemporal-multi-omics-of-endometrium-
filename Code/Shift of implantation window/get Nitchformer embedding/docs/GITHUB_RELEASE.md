# GitHub and manuscript handoff

This directory is prepared as a standalone source-code contribution. Existing
analysis folders and the source data remain separate.

## Files to publish

Publish this directory's scripts, tests, environment files, docs and provenance
records. The supplied `.gitignore` excludes data, results, weights, caches and
large generated binary files. Review staged filenames before committing.

If adding to an existing repository, run from its root:

```bash
git add Script_10_RIF_implantation_Nitchformer_network/
git diff --cached --stat
```

For a new standalone repository, copy this directory's contents into that repo
(including `.gitignore`). The relative paths in README then apply directly.
No repository remote, branch, author identity or license has been assumed.

## Publication metadata to complete

Supply the study's data accession/download URL and any access conditions for
`pro_sec_with_RIF.h5ad`, paper title/DOI, author names and preferred code license.
Those details were not provided, so this package does not invent citation metadata
or grant a license on the authors' behalf. Public availability alone is not an
explicit reuse license. Attribute the official Nicheformer model and its upstream
license separately; model sources/weights are downloaded rather than vendored.

For manuscript archival, tag the tested code version and deposit the release with
a persistent archive/DOI service chosen by the authors. Include input and model
hashes, all executed commands, runtime versions, sampling seeds and any deviations
from the reference workflow with the study results.

The newly added GPU entry point should be run on the target GPU before claiming
full GPU inference reproduction. The full historical preprocessing was not rerun
when this code package was organized; its original reports remain the reference.
