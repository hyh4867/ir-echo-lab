import pytest

from irecho.healpix import ang2pix_nest


@pytest.mark.parametrize("ra,dec,face", [
    (0.0, 0.0, 4), (90.0, 0.0, 5), (180.0, 0.0, 6), (270.0, 0.0, 7),
    (45.0, 60.0, 0), (135.0, 60.0, 1), (225.0, 60.0, 2), (315.0, 60.0, 3),
    (45.0, -60.0, 8), (315.0, -60.0, 11),
])
def test_base_pixels(ra, dec, face):
    """At order 0 the pixel index is the base-face number (Gorski et al. 2005, Fig. 4)."""
    assert ang2pix_nest(0, ra, dec) == face


def test_nested_hierarchy():
    """A pixel's parent at order k-1 is its index divided by four."""
    for ra, dec in [(174.927, 31.909), (25.99, 2.35), (54.19, -36.0), (220.6, -17.25), (10.0, 85.0)]:
        for order in range(1, 8):
            assert ang2pix_nest(order, ra, dec) >> 2 == ang2pix_nest(order - 1, ra, dec)


def test_ngc3786_partition():
    """The unTimely partition healpix_k5=7148 spans RA 172.97-175.78, Dec 30.00-32.79
    (read from the parquet statistics on 2026-11-12) and contains NGC 3786."""
    assert ang2pix_nest(5, 174.927129, 31.9094276) == 7148
    assert ang2pix_nest(5, 174.927129, 31.9094276) >> 10 == 6
