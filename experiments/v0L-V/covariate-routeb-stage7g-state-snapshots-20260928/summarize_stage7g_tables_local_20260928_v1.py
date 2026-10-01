"""Read-only Stage 7G posthoc tables from the hash-verified result archive."""

import csv
from pathlib import Path


BASE = Path(__file__).resolve().parent
ROOT = BASE / "stage7g_20260928_v1"
ANALYSIS = ROOT / "analysis_v1"
OUT = BASE / "posthoc_v1"
FIT_IDS = (
    "covrb7c_01__C1__start06",
    "covrb7c_02__C1__start11",
    "covrb7c_03__C1__start06",
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def number(row, name):
    value = row[name]
    return float(value) if value != "" else float("nan")


def first_sweep(rows, predicate, minimum=0, maximum=499):
    for row in rows:
        sweep = int(row["iteration"])
        if minimum <= sweep <= maximum and predicate(row):
            return sweep
    return ""


def write_csv(path, rows):
    require(not path.exists(), "Refusing overwrite: {}".format(path))
    require(bool(rows), "Empty output: {}".format(path))
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    require((ROOT / "endpoint_parity_v1" / "PARITY_PASS.txt").is_file(),
            "Stage 7G parity marker missing")
    require((ANALYSIS / "ANALYSIS_COMPLETE.txt").is_file(),
            "Stage 7G analysis marker missing")
    require(not OUT.exists(), "Refusing posthoc overwrite")
    state = read_csv(ANALYSIS / "STATE_METRICS_ALL_SWEEPS.csv")
    truth = read_csv(ANALYSIS / "TRUTH_ALIGNMENT_ALL_SWEEPS.csv")
    require(len(state) == len(truth) == 6000, "Incomplete all-sweep tables")
    key = lambda row: (row["fit_id"], row["iteration"], row["block"],
                       row["study"])
    state_by_key = {key(row): row for row in state}
    truth_by_key = {key(row): row for row in truth}
    require(len(state_by_key) == len(truth_by_key) == 6000 and
            set(state_by_key) == set(truth_by_key), "State/truth keys differ")

    events = []
    snapshots = []
    cross_block = []
    selected_sweeps = {0, 1, 2, 5, 10, 20, 25, 30, 40, 50, 55,
                       60, 65, 70, 80, 90, 99, 100, 110, 120,
                       150, 200, 400, 499}
    cross_sweeps = {25, 50, 60, 70, 90, 100, 150, 499}
    for fit_id in FIT_IDS:
        b2 = []
        for iteration in range(500):
            index = (fit_id, str(iteration), "specific", "2")
            x = state_by_key[index]
            y = truth_by_key[index]
            row = dict(x)
            row.update({name: value for name, value in y.items()
                        if name not in ("fit_id", "data_id", "iteration",
                                        "block", "study")})
            b2.append(row)
            if iteration in selected_sweeps:
                keep = (
                    "fit_id", "iteration", "temperature", "ppi_mean",
                    "ppi_max", "support_count_at_0_5", "loading_mean_norm",
                    "loading_second_norm", "score_mean_rms",
                    "score_second_rms", "function_mean_l2",
                    "function_second_l2",
                    "mean_signal_l2_rms_over_variables",
                    "raw_loading_abs_cosine_to_target",
                    "posterior_loading_direction_energy_fraction",
                    "loading_norm_ratio_to_truth",
                    "function_mean_l2_ratio_to_truth",
                    "score_mean_rms_ratio_to_truth",
                    "mean_signal_ratio_to_truth",
                    "loading_direction_interpretable_at_1pct_truth_norm",
                )
                snapshots.append({name: row[name] for name in keep})

        active = lambda row: (
            row["loading_direction_interpretable_at_1pct_truth_norm"] ==
            "TRUE")
        event = {
            "fit_id": fit_id,
            "first_load_norm_below_1pct_truth": first_sweep(
                b2, lambda r: number(r, "loading_norm_ratio_to_truth")
                < 0.01, minimum=1),
            "first_score_mean_below_1pct_truth": first_sweep(
                b2, lambda r: number(r, "score_mean_rms_ratio_to_truth")
                < 0.01, minimum=1),
            "first_ppi_mean_below_0_01": first_sweep(
                b2, lambda r: number(r, "ppi_mean") < 0.01,
                minimum=1),
            "first_mean_signal_below_0_01_truth": first_sweep(
                b2, lambda r: number(r, "mean_signal_ratio_to_truth")
                < 0.01, minimum=1),
            "first_load_norm_above_1_after_25": first_sweep(
                b2, lambda r: number(r, "loading_mean_norm") >= 1,
                minimum=25),
            "first_score_mean_above_0_01_after_25": first_sweep(
                b2, lambda r: number(r, "score_mean_rms") >= 0.01,
                minimum=25),
            "first_mean_signal_above_0_01_truth_after_25": first_sweep(
                b2, lambda r: number(r, "mean_signal_ratio_to_truth")
                >= 0.01, minimum=25),
            "first_mean_signal_above_0_1_truth_after_25": first_sweep(
                b2, lambda r: number(r, "mean_signal_ratio_to_truth")
                >= 0.1, minimum=25),
            "first_cosine_above_0_5_with_active_load": first_sweep(
                b2, lambda r: active(r) and
                number(r, "raw_loading_abs_cosine_to_target") >= 0.5),
            "first_cosine_above_0_9_with_active_load": first_sweep(
                b2, lambda r: active(r) and
                number(r, "raw_loading_abs_cosine_to_target") >= 0.9),
            "first_cosine_above_0_99_with_active_load": first_sweep(
                b2, lambda r: active(r) and
                number(r, "raw_loading_abs_cosine_to_target") >= 0.99),
            "first_ppi_max_above_0_5_after_25": first_sweep(
                b2, lambda r: number(r, "ppi_max") >= 0.5,
                minimum=25),
            "initial_cosine_to_true_B2":
                b2[0]["raw_loading_abs_cosine_to_target"],
            "cosine_at_sweep_60":
                b2[60]["raw_loading_abs_cosine_to_target"],
            "load_norm_at_sweep_60": b2[60]["loading_mean_norm"],
            "signal_ratio_at_sweep_60":
                b2[60]["mean_signal_ratio_to_truth"],
            "endpoint_cosine":
                b2[499]["raw_loading_abs_cosine_to_target"],
            "endpoint_signal_ratio":
                b2[499]["mean_signal_ratio_to_truth"],
        }
        events.append(event)
        for iteration in sorted(cross_sweeps):
            for block, study, label in (("shared", "2", "shared"),
                                        ("specific", "1", "B1")):
                index = (fit_id, str(iteration), block, study)
                x = state_by_key[index]
                y = truth_by_key[index]
                cross_block.append({
                    "fit_id": fit_id,
                    "iteration": iteration,
                    "fitted_block": label,
                    "cosine_to_true_B2":
                        y["raw_loading_abs_cosine_to_true_specific_2"],
                    "loading_mean_norm": x["loading_mean_norm"],
                    "mean_signal_l2_rms_over_variables":
                        x["mean_signal_l2_rms_over_variables"],
                })
    OUT.mkdir()
    write_csv(OUT / "B2_EVENTS_ALL_SWEEPS.csv", events)
    write_csv(OUT / "B2_SCALE_DIRECTION_KEY.csv", snapshots)
    write_csv(OUT / "OTHER_BLOCK_ALIGNMENT_TO_TRUE_B2.csv", cross_block)
    (OUT / "POSTHOC_COMPLETE.txt").write_text(
        "STAGE7G_POSTHOC_COMPLETE\n"
        "read_only_from_verified_stage7g_tables=TRUE\n"
        "new_fit_started=FALSE\n", encoding="utf-8")
    for row in events:
        print(row)


if __name__ == "__main__":
    main()
