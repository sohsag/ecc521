import pytest
from montgomery25519 import *
from utils import AffinePoint, affine_to_proj, proj_to_affine


@pytest.fixture
def G_in_montgomery():
    Gx = K(0x09)
    Gy = K(0x20AE19A1B8A086B4E01EDD2C7748D14C923D4D7E6D7C61B229E9C5A27ECED3D9)
    G_in_montgomery = AffinePoint(Gx, Gy)
    return G_in_montgomery


@pytest.fixture
def G_in_weierstrass(G_in_montgomery):
    return E(*to_weierstrass(G_in_montgomery))


@pytest.fixture
def G_projective_xz(G_in_montgomery) -> ProjectivePointXZ:
    return normal_proj_to_xz(affine_to_proj(G_in_montgomery, K))


def test_xDBL(G_projective_xz, G_in_weierstrass):
    G2 = xDBL(G_projective_xz)
    G2x = xz_to_x_coordinate(G2)
    G2x = x_in_weierstrass(G2x)

    assert (G_in_weierstrass * 2)[0] == G2x


def test_xADD(G_projective_xz, G_in_weierstrass):
    G2 = xDBL(G_projective_xz)
    diff_G_G2 = G_projective_xz  # we know the difference between these are G

    G3 = xADD(G_projective_xz, G2, diff_G_G2)
    G3x = xz_to_x_coordinate(G3)
    G3x = x_in_weierstrass(G3x)

    assert (G_in_weierstrass * 2 + G_in_weierstrass)[0] == G3x


def test_uniform_ladder1(G_projective_xz, G_in_weierstrass):
    G212, _ = constant_time_montgomery_ladder(G_projective_xz, 212)
    G212x = xz_to_x_coordinate(G212)
    G212x = x_in_weierstrass(G212x)

    assert (G_in_weierstrass * 212)[0] == G212x


def test_uniform_ladder2(G_projective_xz, G_in_weierstrass):
    G313, _ = constant_time_montgomery_ladder(G_projective_xz, 313)
    G313x = xz_to_x_coordinate(G313)
    G313x = x_in_weierstrass(G313x)

    assert (G_in_weierstrass * 313)[0] == G313x


def test_y_recovery(G_projective_xz, G_in_weierstrass, G_in_montgomery):
    G313, G314 = constant_time_montgomery_ladder(G_projective_xz, 313)
    G313_full = y_recovery(G_in_montgomery, G313, G314)
    G313_full = proj_to_affine(G313_full)
    G313_full = E(*to_weierstrass(G313_full))
    y = G313_full.xy()[1]
    assert (G_in_weierstrass * 313)[1] == y


def test_montgomery_ladder(G_in_weierstrass, G_projective_xz):
    G313, _ = montgomery_ladder(G_projective_xz, 313)
    G313x = xz_to_x_coordinate(G313)
    G313x = x_in_weierstrass(G313x)

    assert (G_in_weierstrass * 313)[0] == G313x
