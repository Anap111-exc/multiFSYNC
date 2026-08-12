# v0L-V frozen1 binding

This is the first frozen binding of the candidate2 v0L-V design. It binds the
unchanged rc1 snapshot to the accepted Alibaba Cloud Linux R environment.

- Parent candidate: `protocol_manifests_20260812_v3_candidate2`.
- Runtime source commit: `b49448f851e5401cbaf681577ecee2b4bdc72514`.
- Frozen tag: `v0.3.0-v0L-V-frozen1`.
- Data/fit/audit rows: 10/210/72.
- Runtime source hashes: 50/50 verified on the target Linux server.
- Post-security-update Linux smoke: 12/12 PASS, exit 0.
- Status: frozen, but formal execution is not authorized.

The files named `*CANDIDATE2.csv` retain their original names and identifiers
because their scientific contents are copied byte-for-byte from candidate2.
Only the binding status, Linux environment, installed-function signatures, and
freeze lineage are new.

Static verification from the repository root:

```sh
Rscript experiments/v0L-V/frozen1/verify_frozen1_binding.R
```

On the accepted server, with the package installed from `rc1_snapshot`:

```sh
export R_LIBS="$PWD/experiments/v0L-V/rc1_snapshot/_r_test_lib"
Rscript experiments/v0L-V/frozen1/verify_frozen1_binding.R --require-runtime
```

This binding must remain unauthorized. A later explicit start instruction may
create a separate, versioned authorized execution binding; it must not modify
this directory or the frozen tag.
