# Stage 6G local protocol validation

Date: 2026-09-11  
Scope: protocol/configuration consistency only; no data generation, fit,
continuation, scoring or truth access

Command:

```sh
Rscript validate_stage6g_protocol_20260911_v1.R
```

Windows R 4.5.1 could not preserve the Chinese parent path after
`normalizePath()`, so the identical files were checked through a temporary
ASCII `subst` path. The validator itself now deliberately preserves the path
supplied by `--file`, which works both with the temporary Windows mapping and
normal Linux absolute paths.

Result: 17/17 registered checks passed.

- scenarios: 3;
- common candidate configurations: 5;
- holdout splits: 2;
- paired fitting starts: 4;
- implied endpoints: 120;
- configuration/split winners: 30;
- planning CPU time: 96.72 single-core hours;
- all seeds and IDs internally unique;
- all three generating truths are present in the common grid;
- all truth dimensions are inside the currently validated exact upstream
  Jaoua/Splinets template domain;
- complexity proxies reproduce the registered formula exactly;
- execution and formal-paper flags remain false;
- all protocol R files parse.

No ECS connection was made and no package source was changed by this
validation.

