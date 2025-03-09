import pytest
from edwards import *

@pytest.fixture
def G_in_edwards():
    Gx = K(
        0x752cb45c48648b189df90cb2296b2878a3bfd9f42fc6c818ec8bf3c9c0c6203913f6ecc5ccc72434b1ae949d568fc99c6059d0fb13364838aa302a940a2f19ba6c)
    Gy = K(0x0c)
    G_edwards = AffinePoint(Gx, Gy)
    return G_edwards

@pytest.fixture
def G_in_weierstrass(G_in_edwards):
    return E(*to_weierstrass(d, G_in_edwards))

def test_simple_translation(G_in_edwards):
    G = E(*to_weierstrass(d, G_in_edwards))
    point = to_edwards(d, AffinePoint(G[0], G[1]))

    assert point.x == G_in_edwards[0]
    assert point.y == G_in_edwards[1]


def test_addition(G_in_edwards, G_in_weierstrass):
    G_x, G_y = G_in_edwards
    G_proj = affine_to_proj(G_x, G_y)
    TwoG = edwards_addition(G_proj, G_proj)
    affine_edwards_two_g = proj_to_affine(*TwoG)
    affine_weierstrass_two_g = to_weierstrass(d, affine_edwards_two_g)
    assert (G_in_weierstrass * 2)[0] == affine_weierstrass_two_g.x
    assert (G_in_weierstrass * 2)[1] == affine_weierstrass_two_g.y


