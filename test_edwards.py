import pytest
from edwards521 import *


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

@pytest.fixture
def G_projective(G_in_edwards):
    return affine_to_proj(G_in_edwards)


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


def test_doubling(G_projective, G_in_weierstrass):
    TwoG = edwards_doubling(G_projective)
    affine_edwards_two_g = proj_to_affine(TwoG)
    affine_weierstrass_two_g = to_weierstrass(affine_edwards_two_g)

    assert (G_in_weierstrass * 2)[0] == affine_weierstrass_two_g.x
    assert (G_in_weierstrass * 2)[1] == affine_weierstrass_two_g.y


def test_scalar_multiplication(G_projective, G_in_weierstrass):
    Q = edwards_scalar_multiplication(G_projective, 2)

    affine_edwards_q = proj_to_affine(Q)
    affine_weierstrass_q = to_weierstrass(affine_edwards_q)

    assert (G_in_weierstrass * 2)[0] == affine_weierstrass_q.x
    assert (G_in_weierstrass * 2)[1] == affine_weierstrass_q.y


def test_add_inversions(G_projective):
    G_proj_inv = edwards_negation(G_projective)
    neutral_element = edwards_addition(G_projective, G_proj_inv)
    assert neutral_element.X == 0
    assert neutral_element.Y == neutral_element.Z

    neutral_affine = proj_to_affine(neutral_element)
    assert neutral_affine == Edwards_Neutral_Element_Affine  # neutral element for edwards curve in affine coordinates
