"""Unit tests for aggregate_frequency_tables.py — main() and parse_args()."""
import pytest
import pandas as pd

import aggregate_frequency_tables as agg


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

FREQ_COLS = ["CN", "counts_obs", "freq_obs", "freq_theoretical",
             "counts_corrected", "freq_corrected", "detection_limit"]


def _write_freq_tsv(path, rows):
    pd.DataFrame(rows).to_csv(path, index=False, sep="\t")
    return path


def _make_sample_table(path, samples, paths):
    pd.DataFrame({"sample": samples, "path": paths}).to_csv(path, index=False, sep="\t")
    return path


@pytest.fixture()
def two_samples(tmp_path):
    rows_s1 = [dict(zip(FREQ_COLS, [0, 90, 0.9, 1.0, 90.0, 0.45, 0.01]),),
               dict(zip(FREQ_COLS, [1, 10, 0.1, 0.5, 20.0, 0.55, 0.02]))]
    rows_s2 = [dict(zip(FREQ_COLS, [0, 80, 0.8, 1.0, 80.0, 0.40, 0.01])),
               dict(zip(FREQ_COLS, [1, 15, 0.15, 0.5, 30.0, 0.50, 0.02])),
               dict(zip(FREQ_COLS, [2,  5, 0.05, 0.1,  50.0, 0.10, 0.10]))]

    f1 = _write_freq_tsv(tmp_path / "s1.tsv", rows_s1)
    f2 = _write_freq_tsv(tmp_path / "s2.tsv", rows_s2)
    st = _make_sample_table(tmp_path / "samples.tsv", ["sampleA", "sampleB"],
                            [str(f1), str(f2)])
    return dict(
        sample_table=str(st),
        input_files=[str(f1), str(f2)],
        output_tsv=str(tmp_path / "out.tsv"),
        output_xlsx=str(tmp_path / "out.xlsx"),
        rows_s1=rows_s1,
        rows_s2=rows_s2,
    )


@pytest.fixture()
def one_sample(tmp_path):
    rows = [dict(zip(FREQ_COLS, [0, 100, 1.0, 1.0, 100.0, 1.0, 0.01]))]
    f = _write_freq_tsv(tmp_path / "s1.tsv", rows)
    st = _make_sample_table(tmp_path / "samples.tsv", ["sampleX"], [str(f)])
    return dict(
        sample_table=str(st),
        input_files=[str(f)],
        output_tsv=str(tmp_path / "out.tsv"),
        output_xlsx=str(tmp_path / "out.xlsx"),
        rows=rows,
    )


# ---------------------------------------------------------------------------
# main() — row counts and sample labelling
# ---------------------------------------------------------------------------

def test_main_total_row_count(two_samples):
    agg.main(**{k: two_samples[k] for k in ("sample_table", "input_files", "output_tsv", "output_xlsx")})
    result = pd.read_csv(two_samples["output_tsv"], sep="\t")
    assert len(result) == len(two_samples["rows_s1"]) + len(two_samples["rows_s2"])


def test_main_sample_column_values(two_samples):
    agg.main(**{k: two_samples[k] for k in ("sample_table", "input_files", "output_tsv", "output_xlsx")})
    result = pd.read_csv(two_samples["output_tsv"], sep="\t")
    assert set(result["sample"]) == {"sampleA", "sampleB"}
    assert len(result[result["sample"] == "sampleA"]) == len(two_samples["rows_s1"])
    assert len(result[result["sample"] == "sampleB"]) == len(two_samples["rows_s2"])


def test_main_preserves_original_columns(two_samples):
    agg.main(**{k: two_samples[k] for k in ("sample_table", "input_files", "output_tsv", "output_xlsx")})
    result = pd.read_csv(two_samples["output_tsv"], sep="\t")
    for col in FREQ_COLS:
        assert col in result.columns


def test_main_single_sample(one_sample):
    agg.main(**{k: one_sample[k] for k in ("sample_table", "input_files", "output_tsv", "output_xlsx")})
    result = pd.read_csv(one_sample["output_tsv"], sep="\t")
    assert len(result) == len(one_sample["rows"])
    assert list(result["sample"]) == ["sampleX"]


def test_main_numeric_looking_sample_name_kept_as_string(tmp_path):
    """Sample names that look like numbers must not be cast to int."""
    rows = [{"CN": 0, "freq_obs": 1.0}]
    f = _write_freq_tsv(tmp_path / "s.tsv", rows)
    st = _make_sample_table(tmp_path / "samples.tsv", ["001"], [str(f)])
    out_tsv = str(tmp_path / "out.tsv")
    agg.main(str(st), [str(f)], out_tsv, str(tmp_path / "out.xlsx"))
    result = pd.read_csv(out_tsv, sep="\t", dtype={"sample": str})
    assert result["sample"].iloc[0] == "001"


# ---------------------------------------------------------------------------
# main() — TSV output format
# ---------------------------------------------------------------------------

def test_main_tsv_is_tab_separated(two_samples):
    agg.main(**{k: two_samples[k] for k in ("sample_table", "input_files", "output_tsv", "output_xlsx")})
    with open(two_samples["output_tsv"]) as f:
        header = f.readline()
    assert "\t" in header
    assert "," not in header


def test_main_tsv_has_no_index_column(two_samples):
    agg.main(**{k: two_samples[k] for k in ("sample_table", "input_files", "output_tsv", "output_xlsx")})
    result = pd.read_csv(two_samples["output_tsv"], sep="\t")
    assert "Unnamed: 0" not in result.columns


def test_main_creates_tsv_file(two_samples):
    import os
    agg.main(**{k: two_samples[k] for k in ("sample_table", "input_files", "output_tsv", "output_xlsx")})
    assert os.path.isfile(two_samples["output_tsv"])


# ---------------------------------------------------------------------------
# main() — XLSX output
# ---------------------------------------------------------------------------

def test_main_creates_xlsx_file(two_samples):
    import os
    agg.main(**{k: two_samples[k] for k in ("sample_table", "input_files", "output_tsv", "output_xlsx")})
    assert os.path.isfile(two_samples["output_xlsx"])


def test_main_xlsx_sheet_name(two_samples):
    import openpyxl
    agg.main(**{k: two_samples[k] for k in ("sample_table", "input_files", "output_tsv", "output_xlsx")})
    wb = openpyxl.load_workbook(two_samples["output_xlsx"])
    assert "Frequencies" in wb.sheetnames


def test_main_xlsx_matches_tsv(two_samples):
    agg.main(**{k: two_samples[k] for k in ("sample_table", "input_files", "output_tsv", "output_xlsx")})
    tsv = pd.read_csv(two_samples["output_tsv"], sep="\t")
    xlsx = pd.read_excel(two_samples["output_xlsx"], sheet_name="Frequencies")
    assert list(tsv.columns) == list(xlsx.columns)
    assert len(tsv) == len(xlsx)


# ---------------------------------------------------------------------------
# parse_args()
# ---------------------------------------------------------------------------

def test_parse_args_all_arguments():
    args = agg.parse_args([
        "--sample-table", "samples.tsv",
        "--input", "a.tsv", "b.tsv",
        "--output-tsv", "out.tsv",
        "--output-xlsx", "out.xlsx",
    ])
    assert args.sample_table == "samples.tsv"
    assert args.input_files == ["a.tsv", "b.tsv"]
    assert args.output_tsv == "out.tsv"
    assert args.output_xlsx == "out.xlsx"


def test_parse_args_single_input_file():
    args = agg.parse_args([
        "--sample-table", "s.tsv",
        "--input", "only.tsv",
        "--output-tsv", "out.tsv",
        "--output-xlsx", "out.xlsx",
    ])
    assert args.input_files == ["only.tsv"]


def test_parse_args_missing_sample_table_exits():
    with pytest.raises(SystemExit):
        agg.parse_args(["--input", "a.tsv", "--output-tsv", "o.tsv", "--output-xlsx", "o.xlsx"])


def test_parse_args_missing_input_exits():
    with pytest.raises(SystemExit):
        agg.parse_args(["--sample-table", "s.tsv", "--output-tsv", "o.tsv", "--output-xlsx", "o.xlsx"])


def test_parse_args_missing_output_tsv_exits():
    with pytest.raises(SystemExit):
        agg.parse_args(["--sample-table", "s.tsv", "--input", "a.tsv", "--output-xlsx", "o.xlsx"])


def test_parse_args_missing_output_xlsx_exits():
    with pytest.raises(SystemExit):
        agg.parse_args(["--sample-table", "s.tsv", "--input", "a.tsv", "--output-tsv", "o.tsv"])
