# Local validation record

Date: 2026-08-30
Scope: nonformal Windows micro smoke; no pilot-scale data or fit
Formal paper Monte Carlo result: `FALSE`

Validation results:

- all four pilot R files parsed successfully;
- manifest QC passed `10/10` with `2` data rows, `24` G12 rows, and `2`
  pooled bayesSYNC rows;
- the source binding matched all `30/30` registered files before this record
  was added; the final binding is regenerated after this record;
- a fresh external temporary output directory completed generation/sealing,
  two micro G12 fits, two micro pooled fits, truth-free selection of two G12
  winners, smoke-only authorized unsealing, four evaluations, and resource
  summarization;
- terminal fits completed `4/4`, objective-eligible endpoints were `4/4`, and
  evaluation errors were `0/4`;
- metric completeness passed `4/4`: each G12 evaluation produced all 15
  required evaluator/Stage-5A tables and each pooled evaluation produced all 5
  required tables;
- the smoke output contained `SMOKE_COMPLETE.txt` with `status=PASS` and stayed
  outside the Git repository;
- no full-size pilot data, 26-fit pilot, truth unseal, continuation, or formal
  paper Monte Carlo was started by this local validation.

The micro smoke checks software wiring and metric presence only. Its numerical
recovery values are not scientific results and are not retained in Git.
