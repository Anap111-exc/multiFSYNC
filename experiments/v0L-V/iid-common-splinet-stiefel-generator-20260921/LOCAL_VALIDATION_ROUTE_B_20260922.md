# Route-B local validation record

Date: 2026-09-22

- Package R parse: PASS, 27 files.
- Targeted `structured_generator` test: PASS.
- Full `testthat`: PASS, no failures; 15 existing max-iteration warnings.
- Route-B no-fit K audit: PASS; minimum eligible `K=8`.
- Route-B generation-only smoke: PASS.
- Formal data generated: FALSE.
- Model fits started: FALSE.
- Continuations started: FALSE.

The Windows Unicode project path still triggers the known `pkgload/desc`
package-root issue. Tests were run from a byte-for-byte ASCII temporary copy;
the source package in the project directory was not rewritten by the test
workaround.

