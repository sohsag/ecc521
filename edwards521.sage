from utils import (
    Edwards_Neutral_Element_Affine,
    Edwards_Neutral_Element_Projective,
    AffinePoint,
    ProjectivePoint
)
# constants from https://neuromancer.sk/std/other/E-521
p = 0x1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
K = GF(p)
# d is used to construct edwards curve.
d = K(
    0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffa4331)
# This is a weierstrass curve, that is birationally equivalent to the edwards curve with the d value above.
E = EllipticCurve(K, (
    K(-1 / 48) * (1 + 14 * d + d ^ 2),  # a
    K(1 / 864) * (1 + d) * (-1 + 34 * d - d ^ 2)  # b
))

# conversions found https://neuromancer.sk/std/other/Ed25519
def to_weierstrass(P: AffinePoint) -> AffinePoint:
    x = (5 + P.y - 5 * d * P.y - d) / (12 - 12 * P.y)
    y = (1 + P.y - d * P.y - d) / (4 * P.x - 4 * P.x * P.y)
    return AffinePoint(x, y)


def to_edwards(P: AffinePoint) -> AffinePoint:
    u, v = P
    y = (5 - 12 * u - d) / (-12 * u - 1 + 5 * d)
    x = (1 + y - d * y - d) / (4 * v - 4 * v * y)
    return AffinePoint(x, y)

def edwards_addition(P1: ProjectivePoint, P2: ProjectivePoint) -> ProjectivePoint:
    # page 21
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
    Y3 = A * G * (D - C)
    Z3 = G * F

    return ProjectivePoint(X3, Y3, Z3)

def edwards_scalar_multiplication(P: ProjectivePoint, scalar: int) -> ProjectivePoint:
    Q = P
    R = Edwards_Neutral_Element_Projective
    while scalar > 0:
        if scalar % 2 == 1:
            R = edwards_addition(R, Q)
        Q = edwards_addition(Q, Q)
        scalar = scalar // 2

    return R


def edwards_negation(P: ProjectivePoint) -> ProjectivePoint:
    X, Y, Z = P
    return ProjectivePoint(-X, Y, Z)


"""
def montgomery_ladder(P: ProjectivePoint, k: int):
    k_bin = bin(k)[2:] 
    k_bin = k_bin[::-1]  
    R0, R1 = P, edwards_doubling(P)

    for i in range(len(k_bin) - 2, -1, -1): 
        if k_bin[i] == "0":
            R0, R1 = edwards_doubling(R0), edwards_addition(R0, R1)
        else:
            R0, R1 = edwards_addition(R0, R1), edwards_doubling(R1)

    return R0
    """