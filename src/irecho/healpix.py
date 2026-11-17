"""HEALPix NESTED pixel index of a sky position (Gorski et al. 2005, ApJ 622, 759).

Only what is needed to find the partition of an IRSA HEALPix-partitioned
parquet catalogue; it avoids a dependency on healpy.
"""

from __future__ import annotations

import math


def _interleave(ix: int, iy: int) -> int:
    """Bits of ix on even positions, bits of iy on odd positions."""
    out = 0
    for bit in range(32):
        out |= ((ix >> bit) & 1) << (2 * bit)
        out |= ((iy >> bit) & 1) << (2 * bit + 1)
    return out


def ang2pix_nest(order: int, ra_deg: float, dec_deg: float) -> int:
    """NESTED pixel index at resolution ``order`` (nside = 2**order)."""
    nside = 1 << order
    z = math.sin(math.radians(dec_deg))
    za = abs(z)
    tt = (math.radians(ra_deg) % (2.0 * math.pi)) / (0.5 * math.pi)  # in [0, 4)
    if za <= 2.0 / 3.0:  # equatorial region
        temp1 = nside * (0.5 + tt)
        temp2 = nside * z * 0.75
        jp = int(temp1 - temp2)  # index of ascending edge line
        jm = int(temp1 + temp2)  # index of descending edge line
        ifp, ifm = jp >> order, jm >> order
        if ifp == ifm:
            face = ifp | 4
        elif ifp < ifm:
            face = ifp
        else:
            face = ifm + 8
        ix = jm & (nside - 1)
        iy = nside - (jp & (nside - 1)) - 1
    else:  # polar caps
        ntt = min(int(tt), 3)
        tp = tt - ntt
        tmp = nside * math.sqrt(3.0 * (1.0 - za))
        jp = min(int(tp * tmp), nside - 1)
        jm = min(int((1.0 - tp) * tmp), nside - 1)
        if z >= 0:
            face, ix, iy = ntt, nside - jm - 1, nside - jp - 1
        else:
            face, ix, iy = ntt + 8, jp, jm
    return face * nside * nside + _interleave(ix, iy)
