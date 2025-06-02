import pytest
from edwards25519 import (
    to_weierstrass,
    E,
    to_twistededwards,
    K,
    edwards_addition,
    edwards_scalar_multiplication,
    encode_point,
    decode_point,
    generate_keypair,
    sign_message,
    verify_message,
)
from utils import AffinePoint, affine_to_proj, proj_to_affine


@pytest.fixture
def G_in_edwards():
    Gx, Gy = K(0x216936D3CD6E53FEC0A4E231FDD6DC5C692CC7609525A7B2C9562D608F25D51A), K(
        0x6666666666666666666666666666666666666666666666666666666666666658
    )
    G_in_edwards = AffinePoint(Gx, Gy)
    return G_in_edwards


@pytest.fixture
def G_in_weierstrass(G_in_edwards):
    return E(*to_weierstrass(G_in_edwards))


@pytest.fixture
def G_projective(G_in_edwards):
    return affine_to_proj(G_in_edwards, K)


def test_simple_translation(G_in_edwards, G_in_weierstrass):
    p = to_twistededwards(AffinePoint(G_in_weierstrass[0], G_in_weierstrass[1]))

    assert p.x == G_in_edwards[0]
    assert p.y == G_in_edwards[1]


def test_addition(G_projective, G_in_weierstrass):
    TwoG = edwards_addition(G_projective, G_projective)
    affine_edwards_two_g = proj_to_affine(TwoG)
    affine_weierstrass_two_g = to_weierstrass(affine_edwards_two_g)

    assert (G_in_weierstrass * 2)[0] == affine_weierstrass_two_g.x
    assert (G_in_weierstrass * 2)[1] == affine_weierstrass_two_g.y


def test_scalar_multiplication(G_in_edwards, G_in_weierstrass):
    Q = edwards_scalar_multiplication(G_in_edwards, 2)
    affine_weierstrass_q = to_weierstrass(Q)

    assert (G_in_weierstrass * 2)[0] == affine_weierstrass_q.x
    assert (G_in_weierstrass * 2)[1] == affine_weierstrass_q.y


def test_encode_positive(G_in_edwards):
    encoded_point = encode_point(G_in_edwards)
    assert encoded_point[0] & 0x80 == 0


def test_encode_negative(G_in_edwards):
    G_in_edwards = edwards_scalar_multiplication(G_in_edwards, 8)
    encoded_point = encode_point(G_in_edwards)
    assert encoded_point[0] & 0x80 == 128


def test_encode_decode(G_in_edwards):
    encoded_point = encode_point(G_in_edwards)
    decoded_point = decode_point(encoded_point)
    assert G_in_edwards == decoded_point


def test_sign_accept():
    A, k = generate_keypair()
    message = b"test_msg"
    R_encoded_int, S = sign_message(message, k)
    msg = verify_message(message, A, R_encoded_int, S)
    assert msg == "Signature accepted"


def test_sign_decline():
    A, k = generate_keypair()

    message = b"test_msg"
    R_encoded_int, S = sign_message(message, k)

    message2 = b"wrong message"
    with pytest.raises(ValueError):
        verify_message(message2, A, R_encoded_int, S)
