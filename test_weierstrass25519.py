from weierstrass25519 import *
from utils import affine_to_jacob, jacob_to_affine

def test_addition():
    G2 = weierstrass_addition(G, G)
    assert (E(*G) * 2)[0] == G2.x
    assert (E(*G) * 2)[1] == G2.y


def test_scalar_mult():
    G8 = weierstrass_scalar_multiplication(G, 8)
    assert (E(*G) * 8)[0] == G8.x
    assert (E(*G) * 8)[1] == G8.y

def test_fermat_scalar_mult():
    G8 = weierstrass_scalar_multiplication_fermat(G, 8)
    assert (E(*G) * 8)[0] == G8.x
    assert (E(*G) * 8)[1] == G8.y

def test_jacob_addition():
    G_jacob = affine_to_jacob(G, K)
    G2_jacob = weierstrass_jacob_addition(G_jacob, G_jacob)
    G2 = jacob_to_affine(G2_jacob)
    assert (E(*G) * 2)[0] == G2.x
    assert (E(*G) * 2)[1] == G2.y

def test_jacob_scalar_multiplication():
    G_jacob = affine_to_jacob(G, K)
    G8_jacob = weierstrass_jacob_scalar_multiplication(G_jacob, 8)
    G8 = jacob_to_affine(G8_jacob)
    assert (E(*G) * 8)[0] == G8.x
    assert (E(*G) * 8)[1] == G8.y