import json

import numpy as np
import pytest

from irecho import irsa


class FakeResponse:
    def __init__(self, content: bytes):
        self.content = content

    def raise_for_status(self):
        return None


class FakeSession:
    """Stands in for ``requests``; counts calls so cache reuse can be asserted."""

    def __init__(self, content: bytes):
        self.content = content
        self.calls = []

    def get(self, url, params=None, timeout=None, headers=None):
        self.calls.append({"url": url, "params": params})
        return FakeResponse(self.content)


CSV = b"ra,dec,mjd,cc_flags,moon_masked,w1mpro,w1sigmpro\n150.0,20.0,56700.1,0000,00,10.300,0.020\n150.0,20.0,56700.2,0h00,01,10.310,\n"


@pytest.fixture(autouse=True)
def no_throttle(monkeypatch):
    monkeypatch.setattr(irsa, "MIN_SECONDS_BETWEEN_QUERIES", 0.0)


def test_cone_adql_is_deterministic_and_complete():
    adql = irsa.cone_adql("neowiser_p1bs_psd", ("ra", "dec", "mjd"), 174.927129, 31.9094276, 3.0)
    assert adql == (
        "SELECT ra, dec, mjd FROM neowiser_p1bs_psd WHERE CONTAINS(POINT('ICRS', ra, dec), "
        "CIRCLE('ICRS', 174.9271290, 31.9094276, 0.0008333)) = 1 ORDER BY mjd"
    )


def test_read_csv_table_types_and_nulls(tmp_path):
    path = tmp_path / "t.csv"
    path.write_bytes(CSV)
    table = irsa.read_csv_table(path)
    assert table["cc_flags"].tolist() == ["0000", "0h00"]  # kept as text
    assert table["moon_masked"].tolist() == ["00", "01"]  # leading zeros preserved
    assert table["w1mpro"].tolist() == [10.3, 10.31]
    assert np.isnan(table["w1sigmpro"][1])  # empty cell -> NaN


def test_fetch_csv_writes_cache_meta_and_manifest(tmp_path):
    session = FakeSession(CSV)
    path = tmp_path / "neowise" / "x.csv"
    manifest = tmp_path / "manifest.json"
    out = irsa.fetch_csv("SELECT 1", path, session=session, manifest_path=manifest)
    assert out == path and path.read_bytes() == CSV
    assert session.calls[0]["params"] == {"QUERY": "SELECT 1", "FORMAT": "csv"}
    meta = json.loads((tmp_path / "neowise" / "x.csv.meta.json").read_text())
    assert meta["n_rows"] == 2
    assert meta["sha256"] == irsa.sha256_bytes(CSV) == irsa.sha256_file(path)
    assert meta["adql"] == "SELECT 1"
    assert json.loads(manifest.read_text())["x.csv"]["sha256"] == meta["sha256"]


def test_fetch_csv_reuses_cache_without_network(tmp_path):
    session = FakeSession(CSV)
    path = tmp_path / "x.csv"
    irsa.fetch_csv("SELECT 1", path, session=session)
    irsa.fetch_csv("SELECT 1", path, session=session)
    assert len(session.calls) == 1
    irsa.fetch_csv("SELECT 1", path, session=session, refresh=True)
    assert len(session.calls) == 2


def test_fetch_csv_refuses_to_reuse_cache_for_a_different_query(tmp_path):
    session = FakeSession(CSV)
    path = tmp_path / "x.csv"
    irsa.fetch_csv("SELECT 1", path, session=session)
    with pytest.raises(RuntimeError, match="different query"):
        irsa.fetch_csv("SELECT 2", path, session=session)


def test_fetch_csv_rejects_error_documents(tmp_path, monkeypatch):
    monkeypatch.setattr(irsa.time, "sleep", lambda s: None)
    session = FakeSession(b"<VOTABLE><INFO name='QUERY_STATUS' value='ERROR'/></VOTABLE>")
    with pytest.raises(RuntimeError, match="non-CSV"):
        irsa.fetch_csv("SELECT 1", tmp_path / "x.csv", session=session, retries=2)
    assert not (tmp_path / "x.csv").exists()
    assert len(session.calls) == 2


def test_count_csv_rows():
    assert irsa.count_csv_rows(CSV.decode()) == 2
    assert irsa.count_csv_rows("a,b\n") == 0


# --- several positions in one query (TAP upload) -------------------------------
class FakePostSession:
    def __init__(self, content: bytes):
        self.content = content
        self.calls = []

    def post(self, url, data=None, files=None, timeout=None, headers=None):
        self.calls.append({"url": url, "data": data, "files": files})
        return FakeResponse(self.content)


UPLOAD_CSV = (b"pid,ra,dec,mjd,cc_flags,w1mpro\n"
              b"1,150.0,20.0,56700.2,0000,10.31\n"
              b"0,151.0,21.0,56900.1,0000,11.00\n"
              b"1,150.0,20.0,56700.1,0h00,10.30\n")


def test_positions_votable_and_upload_adql():
    vot = irsa.positions_votable([(0, 174.927129, 31.9094276), (7, 10.0, -5.0)])
    assert "<TR><TD>0</TD><TD>174.9271290</TD><TD>31.9094276</TD></TR>" in vot
    assert "<TR><TD>7</TD><TD>10.0000000</TD><TD>-5.0000000</TD></TR>" in vot
    assert vot.count("<FIELD") == 3
    adql = irsa.upload_cone_adql("neowiser_p1bs_psd", ("ra", "mjd"), 3.0, "t.w1mpro < 12")
    assert adql == (
        "SELECT p.pid, t.ra, t.mjd FROM neowiser_p1bs_psd AS t, TAP_UPLOAD.pos AS p WHERE "
        "CONTAINS(POINT('ICRS', t.ra, t.dec), CIRCLE('ICRS', p.pra, p.pdec, 0.0008333)) = 1 "
        "AND t.w1mpro < 12"
    )


def test_sort_csv_rows_is_numeric_and_stable_in_content():
    out = irsa.sort_csv_rows(UPLOAD_CSV, ("pid", "mjd")).decode().splitlines()
    assert out[0] == "pid,ra,dec,mjd,cc_flags,w1mpro"
    assert [line.split(",")[0] for line in out[1:]] == ["0", "1", "1"]
    assert [line.split(",")[3] for line in out[1:]] == ["56900.1", "56700.1", "56700.2"]
    # numeric, not lexical: 9 sorts before 10
    assert irsa.sort_csv_rows(b"pid,x\n10,a\n9,b\n", ("pid",)).decode().splitlines()[1] == "9,b"


def test_fetch_upload_csv_caches_sorts_and_records(tmp_path):
    session = FakePostSession(UPLOAD_CSV)
    path = tmp_path / "field.csv"
    manifest = tmp_path / "manifest.json"
    positions = [(0, 151.0, 21.0), (1, 150.0, 20.0)]
    irsa.fetch_upload_csv("tbl", ("ra", "dec", "mjd", "cc_flags", "w1mpro"), positions, 3.0,
                          path, sort_keys=("pid", "mjd"), manifest_path=manifest, session=session)
    assert len(session.calls) == 1
    assert session.calls[0]["data"]["UPLOAD"] == "pos,param:pos"
    meta = json.loads((tmp_path / "field.csv.meta.json").read_text())
    assert meta["n_rows"] == 3 and meta["upload"] == [[0, 151.0, 21.0], [1, 150.0, 20.0]]
    assert meta["sha256"] == irsa.sha256_file(path) and meta["sorted_by"] == ["pid", "mjd"]
    assert json.loads(manifest.read_text())["field.csv"]["sha256"] == meta["sha256"]
    # second call: no request; different positions: refused
    irsa.fetch_upload_csv("tbl", ("ra", "dec", "mjd", "cc_flags", "w1mpro"), positions, 3.0,
                          path, sort_keys=("pid", "mjd"), session=session)
    assert len(session.calls) == 1
    with pytest.raises(RuntimeError, match="different query or position list"):
        irsa.fetch_upload_csv("tbl", ("ra", "dec", "mjd", "cc_flags", "w1mpro"),
                              [(0, 151.0, 21.0)], 3.0, path, session=session)
    parts = irsa.split_by_pid(irsa.read_csv_table(path))
    assert sorted(parts) == [0, 1]
    assert parts[1]["mjd"].tolist() == [56700.1, 56700.2] and "pid" not in parts[1]
    assert parts[0]["w1mpro"].tolist() == [11.0]


def test_fetch_upload_csv_rejects_non_csv(tmp_path):
    session = FakePostSession(b"<VOTABLE><INFO name='QUERY_STATUS' value='ERROR'/></VOTABLE>")
    with pytest.raises(RuntimeError, match="non-CSV"):
        irsa.fetch_upload_csv("tbl", ("ra",), [(0, 1.0, 2.0)], 3.0, tmp_path / "x.csv",
                              session=session, retries=1)
    assert not (tmp_path / "x.csv").exists()
