import gzip
from pathlib import Path

from analysis.tier3_common import _open_text


def test_open_text_detects_gzip_magic_without_filename_suffix(tmp_path: Path):
    object_path = tmp_path / ("a" * 64)
    object_path.write_bytes(gzip.compress(b"##gff-version 3\n", mtime=0))

    with _open_text(object_path) as handle:
        assert handle.read() == "##gff-version 3\n"


def test_open_text_keeps_suffixless_plain_text(tmp_path: Path):
    object_path = tmp_path / ("b" * 64)
    object_path.write_text("plain\n")

    with _open_text(object_path) as handle:
        assert handle.read() == "plain\n"
