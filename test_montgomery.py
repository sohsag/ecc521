import pytest
from montgomery25519 import *
from utils import (
    AffinePoint,
    Montgomery_Neutral_Element_Projective,
    affine_to_proj,
    proj_to_affine
)


@pytest.fixture
def G_in_montgomery():
    Gx = K(0x09)
    Gy = K(0x20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9)
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
    diff_G_G2 = G_projective_xz # we know the difference between these are G

    G3 = xADD(G_projective_xz, G2, diff_G_G2) 
    G3x = xz_to_x_coordinate(G3)
    G3x = x_in_weierstrass(G3x)

    assert (G_in_weierstrass * 2 + G_in_weierstrass)[0] == G3x

def test_uniform_ladder(G_projective_xz, G_in_weierstrass):
    G8 = ladder(G_projective_xz, 8)
    G8x = xz_to_x_coordinate(G8)
    G8x = x_in_weierstrass(G8x)

    assert (G_in_weierstrass * 8)[0] == G8x



