import pytest
from edwards import *


@pytest.fixture
def G_in_edwards():
    Gx = K(
        0x752cb45c48648b189df90cb2296b2878a3bfd9f42fc6c818ec8bf3c9c0c6203913f6ecc5ccc72434b1ae949d568fc99c6059d0fb13364838aa302a940a2f19ba6c)
    Gy = K(0x0c)
    G_in_edwards = AffinePoint(Gx, Gy)
    return G_in_edwards


@pytest.fixture
def G_in_weierstrass(G_in_edwards):
    return E(*to_weierstrass(G_in_edwards))


def test_simple_translation(G_in_edwards):
    G = E(*to_weierstrass(G_in_edwards))
    point = to_edwards(AffinePoint(G[0], G[1]))

    assert point.x == G_in_edwards[0]
    assert point.y == G_in_edwards[1]


def test_addition(G_in_edwards, G_in_weierstrass):
    G_proj = affine_to_proj(G_in_edwards)
    TwoG = edwards_addition(G_proj, G_proj)
    affine_edwards_two_g = proj_to_affine(TwoG)
    affine_weierstrass_two_g = to_weierstrass(affine_edwards_two_g)

    assert (G_in_weierstrass * 2)[0] == affine_weierstrass_two_g.x
    assert (G_in_weierstrass * 2)[1] == affine_weierstrass_two_g.y


def test_doubling(G_in_edwards, G_in_weierstrass):
    G_proj = affine_to_proj(G_in_edwards)
    TwoG = edwards_doubling(G_proj)
    affine_edwards_two_g = proj_to_affine(TwoG)
    affine_weierstrass_two_g = to_weierstrass(affine_edwards_two_g)

    assert (G_in_weierstrass * 2)[0] == affine_weierstrass_two_g.x
    assert (G_in_weierstrass * 2)[1] == affine_weierstrass_two_g.y


def test_scalar_multiplication(G_in_edwards):
    G = E(*to_weierstrass(G_in_edwards))
    G_proj = affine_to_proj(G_in_edwards)
    Q = edwards_scalar_multiplication(G_proj, 2)
    print(Q)
    affine_edwards_q = proj_to_affine(Q)
    affine_weierstrass_q = to_weierstrass(affine_edwards_q)

    assert (G * 2)[0] == affine_weierstrass_q.x
    assert (G * 2)[1] == affine_weierstrass_q.y


def test_add_inversions(G_in_edwards):
    G = E(*to_weierstrass(G_in_edwards))
    G_proj = affine_to_proj(G_in_edwards)
    G_proj_inv = edwards_inversion(G_proj)
    neutral_element = edwards_addition(G_proj, G_proj_inv)
    assert neutral_element.X == 0
    assert neutral_element.Y == neutral_element.Z

    neutral_affine = proj_to_affine(neutral_element)
    assert neutral_affine == Edwards_Neutral_Element_Affine  # neutral element for edwards curve in affine coordinates
