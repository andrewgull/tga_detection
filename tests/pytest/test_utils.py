"""Unit tests for utils.py — get_sample_path() and tsv2dict()."""
import types

import pandas as pd
import pytest

import utils


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _wc(sample):
    """Minimal wildcard-like object with a .sample attribute."""
    return types.SimpleNamespace(sample=sample)


def _make_df(samples, paths):
    return pd.DataFrame({"sample": samples, "path": paths})


# ---------------------------------------------------------------------------
# get_sample_path()
# ---------------------------------------------------------------------------

def test_get_sample_path_returns_correct_path():
    df = _make_df(["s1", "s2"], ["/data/s1.fastq", "/data/s2.fastq"])
    assert utils.get_sample_path(_wc("s1"), df) == "/data/s1.fastq"
    assert utils.get_sample_path(_wc("s2"), df) == "/data/s2.fastq"


def test_get_sample_path_missing_sample_raises_value_error():
    df = _make_df(["s1"], ["/data/s1.fastq"])
    with pytest.raises(ValueError, match="s99"):
        utils.get_sample_path(_wc("s99"), df)


def test_get_sample_path_empty_dataframe_raises_value_error():
    df = _make_df([], [])
    with pytest.raises(ValueError):
        utils.get_sample_path(_wc("s1"), df)


def test_get_sample_path_returns_string():
    df = _make_df(["s1"], ["/data/s1.fastq"])
    result = utils.get_sample_path(_wc("s1"), df)
    assert isinstance(result, str)


def test_get_sample_path_duplicate_sample_returns_first():
    df = _make_df(["s1", "s1"], ["/first.fastq", "/second.fastq"])
    assert utils.get_sample_path(_wc("s1"), df) == "/first.fastq"


def test_get_sample_path_numeric_looking_sample_name():
    df = _make_df(["001", "002"], ["/a.fastq", "/b.fastq"])
    assert utils.get_sample_path(_wc("001"), df) == "/a.fastq"


# ---------------------------------------------------------------------------
# tsv2dict()
# ---------------------------------------------------------------------------

def test_tsv2dict_basic(tmp_path):
    tsv = tmp_path / "params.tsv"
    tsv.write_text("param\tvalue\nalpha\t1\nbeta\thello\n")
    result = utils.tsv2dict(str(tsv))
    assert result == {"alpha": "1", "beta": "hello"}


def test_tsv2dict_returns_dict(tmp_path):
    tsv = tmp_path / "params.tsv"
    tsv.write_text("param\tvalue\nkey\tval\n")
    assert isinstance(utils.tsv2dict(str(tsv)), dict)


def test_tsv2dict_values_are_strings(tmp_path):
    """Numeric-looking values must be kept as strings (dtype enforcement)."""
    tsv = tmp_path / "params.tsv"
    tsv.write_text("param\tvalue\ncount\t42\nrate\t0.5\n")
    result = utils.tsv2dict(str(tsv))
    assert result["count"] == "42"
    assert result["rate"] == "0.5"
    assert all(isinstance(v, str) for v in result.values())


def test_tsv2dict_single_entry(tmp_path):
    tsv = tmp_path / "params.tsv"
    tsv.write_text("param\tvalue\nonly\tone\n")
    assert utils.tsv2dict(str(tsv)) == {"only": "one"}


def test_tsv2dict_empty_file_returns_empty_dict(tmp_path):
    tsv = tmp_path / "params.tsv"
    tsv.write_text("param\tvalue\n")
    assert utils.tsv2dict(str(tsv)) == {}


def test_tsv2dict_param_keys_are_strings(tmp_path):
    """Numeric-looking param names must also be kept as strings."""
    tsv = tmp_path / "params.tsv"
    tsv.write_text("param\tvalue\n1\tabc\n2\txyz\n")
    result = utils.tsv2dict(str(tsv))
    assert "1" in result
    assert "2" in result
