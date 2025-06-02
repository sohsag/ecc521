import pytest
from edwards521 import *
from utils import AffinePoint, affine_to_proj, proj_to_affine


@pytest.fixture
def G_in_edwards():
    Gx = K(
        0x752CB45C48648B189DF90CB2296B2878A3BFD9F42FC6C818EC8BF3C9C0C6203913F6ECC5CCC72434B1AE949D568FC99C6059D0FB13364838AA302A940A2F19BA6C
    )
    Gy = K(0x0C)
    G_in_edwards = AffinePoint(Gx, Gy)
    return G_in_edwards


@pytest.fixture
def G_in_weierstrass(G_in_edwards):
    return E(*to_weierstrass(G_in_edwards))


@pytest.fixture
def G_projective(G_in_edwards):
    return affine_to_proj(G_in_edwards, K)


def test_simple_translation(G_in_edwards, G_in_weierstrass):
    point = to_edwards(AffinePoint(G_in_weierstrass[0], G_in_weierstrass[1]))

    assert point.x == G_in_edwards[0]
    assert point.y == G_in_edwards[1]


def test_addition(G_projective, G_in_weierstrass):
    TwoG = edwards_addition(G_projective, G_projective)
    affine_edwards_two_g = proj_to_affine(TwoG)
    affine_weierstrass_two_g = to_weierstrass(affine_edwards_two_g)

    assert (G_in_weierstrass * 2)[0] == affine_weierstrass_two_g.x
    assert (G_in_weierstrass * 2)[1] == affine_weierstrass_two_g.y


def test_scalar_multiplication(G_in_edwards, G_in_weierstrass):
    TwoG = edwards_scalar_multiplication(G_in_edwards, 8)

    affine_weierstrass_two_g = to_weierstrass(TwoG)

    assert (G_in_weierstrass * 8)[0] == affine_weierstrass_two_g.x
    assert (G_in_weierstrass * 8)[1] == affine_weierstrass_two_g.y
