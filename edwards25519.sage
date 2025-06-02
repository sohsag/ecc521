from utils import (
    Edwards_Neutral_Element_Affine,
    Edwards_Neutral_Element_Projective,
    AffinePoint,
    ProjectivePoint,
    affine_to_proj,
    proj_to_affine,
)
import os
import hashlib
# Constants found https://neuromancer.sk/std/other/Ed25519
p = 0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFED
K = GF(p)
a = K(0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEC)
d = K(0x52036CEE2B6FFE738CC740797779E89800700A4D4141D8AB75EB4DCA135978A3)
E = EllipticCurve(
    K,
    (
        K(-1 / 48) * (a ^ 2 + 14 * a * d + d ^ 2),
        K(1 / 864) * (a + d) * (-a ^ 2 + 34 * a * d - d ^ 2),
    ),
)
E.set_order(0x1000000000000000000000000000000014DEF9DEA2F79CD65812631A5CF5D3ED * 0x08)

# Conversions found https://neuromancer.sk/std/other/Ed25519
def to_weierstrass(P: AffinePoint) -> AffinePoint:
    return AffinePoint(
        (5 * a + a * P.y - 5 * d * P.y - d) / (12 - 12 * P.y),
        (a + a * P.y - d * P.y - d) / (4 * P.x - 4 * P.x * P.y),
    )


def to_twistededwards(P: AffinePoint) -> AffinePoint:
    u, v = P
    y = (5 * a - 12 * u - d) / (-12 * u - a + 5 * d)
    x = (a + a * y - d * y - d) / (4 * v - 4 * v * y)
    return AffinePoint(x, y)


def edwards_addition(P1: ProjectivePoint, P2: ProjectivePoint) -> ProjectivePoint:
    # Youssef El Housni. Edwards curves. 2018. page 21
    X1, Y1, Z1 = P1
    X2, Y2, Z2 = P2

    A = Z1 * Z2
    B = A * A
    C = X1 * X2
    D = Y1 * Y2
    E = d * C * D
    F = B - E
    G = B + E

    X3 = A * F * ((X1 + Y1) * (X2 + Y2) - C - D)
    Y3 = A * G * (D - a * C)
    Z3 = G * F

    return ProjectivePoint(X3, Y3, Z3)


Edwards_Neutral_Element_Projective = ProjectivePoint(K(0), K(1), K(1))


def edwards_scalar_multiplication(P: AffinePoint, scalar: int) -> AffinePoint:
    # Introduction to Mathematical Cryptography table 6.3
    P_proj = affine_to_proj(P, K)
    Q = P_proj
    R = Edwards_Neutral_Element_Projective
    while scalar > 0:
        if scalar % 2 == 1:
            R = edwards_addition(R, Q)
        Q = edwards_addition(Q, Q)
        scalar = scalar // 2

    return proj_to_affine(R)


Gx, Gy = K(0x216936D3CD6E53FEC0A4E231FDD6DC5C692CC7609525A7B2C9562D608F25D51A), K(
    0x6666666666666666666666666666666666666666666666666666666666666658
)
B = AffinePoint(Gx, Gy)
l = 2 ^ 252 + 27742317777372353535851937790883648493  # B's group order


def encode_point(P: AffinePoint) -> bytes:
    # Implemented according to https://datatracker.ietf.org/doc/html/rfc8032#section-5.1.2
    x, y = P

    y_bytes = list(y.to_bytes("little"))
    if x.lift() & 1:
        y_bytes[31] |= 0x80
    return bytes(y_bytes)


def decode_point(encoded_bytes: bytes) -> AffinePoint:
    # Implemented according to https://datatracker.ietf.org/doc/html/rfc8032#section-5.1.3
    encoded = bytearray(encoded_bytes)
    x_0 = (encoded[31] >> 7) & 1
    encoded[31] &= 0x7F
    y_int = int.from_bytes(encoded, "little")

    p = 2 ^ 255 - 19
    if y_int >= p:
        raise ValueError("Decoding fails: y >= p")

    y = K(y_int)

    u = y ^ 2 - 1
    d = -121665 * pow(121666, -1, p)
    v = d * y ^ 2 + 1

    pow_result = pow(u * v ^ 7, (p - 5) // 8, p)
    x = u * v ^ 3 * pow_result

    v_x2 = v * x ^ 2

    if v_x2 == u:
        pass
    elif v_x2 == (-u):
        x = x * pow(2, (p - 1) // 4, p)
    else:
        raise ValueError("Decoding fails: no square root exists")

    if x == 0 and x_0 == 1:
        raise ValueError("Decoding fails: x = 0 and x_0 = 1")

    if x.lift() % 2 != x_0:
        x = p - x

    return AffinePoint(x, y)


def generate_keypair() -> tuple[AffinePoint, bytes]:
    # Implemented accordingly to "High-speed high-security signatures"
    k = os.urandom(32)
    b_least_significant_bits = hashlib.sha512(k).digest()[:32]
    a = int.from_bytes(b_least_significant_bits, "little")

    A = edwards_scalar_multiplication(B, a)
    return A, k


def sign_message(m, k) -> tuple[int, int]:
    # Implemented accordingly to "High-speed high-security signatures"
    secret_key_hash = hashlib.sha512(k).digest()
    b_least_significant_bits = secret_key_hash[:32]

    a = int.from_bytes(b_least_significant_bits, "little") % l

    A = edwards_scalar_multiplication(B, a)
    A_encoded = encode_point(A)
    m_prefix = secret_key_hash[32:]
    r = hashlib.sha512(m_prefix + m).digest()
    r_int = int.from_bytes(r, "little") % l

    R = edwards_scalar_multiplication(B, r_int)
    R_encoded = encode_point(R)

    s = (
        r_int
        + (
            int.from_bytes(hashlib.sha512(R_encoded + A_encoded + m).digest(), "little")
            * a
        )
    ) % l
    return int.from_bytes(R_encoded, "little"), s


def verify_message(m, A, R_encoded_int, s) -> str:
    # Implemented accordingly to "High-speed high-security signatures"
    R_encoded = int.to_bytes(R_encoded_int, 32, "little")
    R = decode_point(R_encoded)
    lhs = edwards_scalar_multiplication(edwards_scalar_multiplication(B, s), 8)
    HAM = (
        int.from_bytes(
            hashlib.sha512(R_encoded + encode_point(A) + m).digest(), "little"
        )
        % l
    )
    rhs_1 = edwards_scalar_multiplication(edwards_scalar_multiplication(A, HAM), 8)
    rhs_2 = edwards_scalar_multiplication(R, 8)
    rhs_1 = affine_to_proj(rhs_1, K)
    rhs_2 = affine_to_proj(rhs_2, K)
    rhs = edwards_addition(rhs_1, rhs_2)
    rhs = proj_to_affine(rhs)
    if lhs != rhs:
        raise ValueError("Signature not valid")

    return "Signature accepted"
