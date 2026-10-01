#!/usr/bin/env python3
"""Read only reduction of the frozen Stage 7I sweeps 51--60 CSV traces.

The script verifies the blind file bindings and never reads sealed truth or
starts a fit. It retains row order because the block trace follows CAVI order.
"""

import csv
import hashlib
import sys
from collections import Counter
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, fieldnames, data):
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: SCRIPT RESULTS_PARENT OUTPUT_DIR")
    parent = Path(sys.argv[1]).resolve(strict=True)
    output = Path(sys.argv[2])
    if output.exists():
        raise SystemExit(f"Refusing existing output: {output}")

    freeze = parent / "stage7i_blind_freeze_20260929_v1"
    marker = (freeze / "BLIND_ENDPOINTS_FROZEN.txt").read_text(encoding="utf-8")
    assert marker.startswith("STAGE7I_BLIND_ENDPOINTS_FROZEN\n")
    assert "truth_bundle_opened=FALSE" in marker
    binding_path = freeze / "BLIND_FILE_HASHES.csv"
    expected_binding_hash = next(
        line.split("=", 1)[1] for line in marker.splitlines()
        if line.startswith("file_hashes_sha256=")
    )
    assert sha256(binding_path) == expected_binding_hash
    bindings = {row["path"]: row["sha256"] for row in rows(binding_path)}

    branch_roots = {
        "U": "stage7i_u_20260929_v1/u_full_v1",
        "A": "stage7i_as_20260929_v1/A_full_v1",
        "S": "stage7i_as_20260929_v1/S_full_v1",
    }
    selected_scale = [
        "stage", "block", "iteration", "temperature", "elbo",
        "objective_scope", "observed_rms", "fitted_rms",
        "mean_contribution_rms", "beta_contribution_rms",
        "shared_contribution_rms", "specific_contribution_rms",
        "expected_rss_per_obs", "mean_residual_sse", "noise_precision_median",
        "mu_rms", "beta_rms", "phi_rms", "psi_rms", "zeta_rms",
        "xi_rms", "a_rms", "b_rms", "xi_score_q_second_moment",
        "specific_slab_loading_q_second_moment",
        "specific_effective_loading_q_second_moment",
        "dense_psi_mean_rms", "dense_psi_qsecond_rms",
        "ppi_mean", "ppi_gt_half", "specific_ppi_mean",
        "elbo_data_likelihood", "elbo_score_xi", "elbo_osullivan_psi",
        "elbo_spike_slab_specific",
    ]
    selected_driver = [
        "sequence", "iteration", "temperature", "annealing", "phase",
        "block", "study", "factor", "component", "data_driver_norm",
        "data_envelope", "data_coherence", "linear_driver_norm",
        "precision_trace", "precision_frobenius", "prior_precision_trace",
        "driver_precision_ratio", "posterior_mean_norm", "slab_mean_norm",
        "slab_variance_mean", "unrestricted_ppi_min",
        "unrestricted_ppi_mean", "unrestricted_ppi_max", "ppi_min",
        "ppi_mean", "ppi_max", "log_prior_odds_mean", "gate_active",
    ]

    binding_rows = []
    inventory = []
    scale_out = []
    driver_out = []
    for branch, base in branch_roots.items():
        for kind in ("SCALE_TRACE.csv", "DRIVER_TRACE.csv"):
            relative = f"{base}/{kind}"
            assert relative in bindings, f"Unbound trace: {relative}"
            path = parent / relative
            actual = sha256(path)
            assert actual == bindings[relative], f"SHA256 mismatch: {relative}"
            binding_rows.append({"branch": branch, "path": relative,
                                 "sha256": actual})
            data = rows(path)
            selected = [row for row in data
                        if row["iteration"].isdigit()
                        and 51 <= int(row["iteration"]) <= 60]
            assert selected, f"Missing 51:60 trace: {relative}"
            present = set(int(row["iteration"]) for row in selected)
            assert present == set(range(51, 61)), f"Incomplete 51:60: {relative}"
            count = Counter(row["block"] for row in selected)
            inventory.extend({"branch": branch, "trace": kind, "block": block,
                              "rows_51_60": n} for block, n in sorted(count.items()))
            if kind == "SCALE_TRACE.csv":
                assert set(selected_scale).issubset(data[0]), "Scale trace schema changed"
                scale_out.extend({"branch": branch, "row_in_trace": i,
                                  **{field: row[field] for field in selected_scale}}
                                 for i, row in enumerate(selected, 1))
            else:
                assert set(selected_driver).issubset(data[0]), "Driver trace schema changed"
                driver_out.extend({"branch": branch, "row_in_trace": i,
                                   **{field: row[field] for field in selected_driver}}
                                  for i, row in enumerate(selected, 1)
                                  if row["study"] == "2" and row["factor"] == "1")

    assert driver_out, "No study-2 factor-1 driver rows found"
    output.mkdir(mode=0o700)
    write_csv(output / "TRACE_INPUT_SHA256.csv",
              ["branch", "path", "sha256"], binding_rows)
    write_csv(output / "TRACE_BLOCK_INVENTORY.csv",
              ["branch", "trace", "block", "rows_51_60"], inventory)
    write_csv(output / "SCALE_BLOCKS_51_60.csv",
              ["branch", "row_in_trace", *selected_scale], scale_out)
    write_csv(output / "DRIVER_STUDY2_FACTOR1_51_60.csv",
              ["branch", "row_in_trace", *selected_driver], driver_out)
    (output / "TRACE_REDUCTION_COMPLETE.txt").write_text(
        "status=COMPLETE\ntruth_opened=FALSE\nnew_fit_started=FALSE\n"
        f"scale_rows={len(scale_out)}\ndriver_b2_rows={len(driver_out)}\n"
        f"script_sha256={sha256(Path(__file__))}\n",
        encoding="utf-8",
    )
    print(f"STAGE7I_TRACE_REDUCTION_COMPLETE scale={len(scale_out)} "
          f"driver_b2={len(driver_out)}")


if __name__ == "__main__":
    main()
