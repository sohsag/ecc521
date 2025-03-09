from sage.all import *
from collections import namedtuple
AffinePoint = namedtuple('AffinePoint', ['x', 'y'])
ProjectivePoint = namedtuple('ProjectivePoint', ['X', 'Y', 'Z'])

# conversions found https://neuromancer.sk/std/other/Ed25519
def to_weierstrass(d, p: AffinePoint):
    x = (5 + p.y - 5 * d * p.y - d) / (12 - 12 * p.y)
    y = (1 + p.y - d * p.y - d) / (4 * p.x - 4 * p.x * p.y)
    return AffinePoint(x, y)


def to_edwards(d, p: AffinePoint):
    u, v = p
    y = (5 - 12 * u - d) / (-12 * u - 1 + 5 * d)
    x = (1 + y - d * y - d) / (4 * v - 4 * v * y)
    return AffinePoint(x, y)


def affine_to_proj(x, y):
    return ProjectivePoint(x, y, K(1))


def proj_to_affine(X, Y, Z):
    return AffinePoint(X / Z, Y / Z)

# constants from https://neuromancer.sk/std/other/E-521
p = 0x1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
K = GF(p)
# d is currently in edwards form
d = K(
    0x1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffa4331)
# This is a weierstrass curve, that is birationally equivalent to the edwards curve with the d value above.
E = EllipticCurve(K, (
    K(-1 / 48) * (1 + 14 * d + d ^ 2),  # a
    K(1 / 864) * (1 + d) * (-1 + 34 * d - d ^ 2)  # b
))


#G = E(*to_weierstrass(d, G_edwards))


def edwards_addition(P1, P2):
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


