from utils import ProjectivePoint, AffinePoint, ProjectivePointXZ
# Constants found https://neuromancer.sk/std/other/Curve25519
p = 0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFED
K = GF(p)
A = K(0x76D06)
B = K(0x01)
E = EllipticCurve(K, ((3 - A ^ 2) / (3 * B ^ 2), (2 * A ^ 3 - 9 * A) / (27 * B ^ 3)))

# Conversions found https://neuromancer.sk/std/other/Curve25519
def to_weierstrass(P: AffinePoint) -> AffinePoint:
    return AffinePoint(P.x / B + A / (3 * B), P.y / B)


def to_montgomery(P: AffinePoint) -> AffinePoint:
    return AffinePoint(B * (P.x - A / (3 * B)), B * P.y)


def normal_proj_to_xz(P: ProjectivePoint) -> ProjectivePointXZ:
    return ProjectivePointXZ(P.X, P.Z)


def x_in_weierstrass(x: int) -> int:
    return x / B + A / (3 * B)


def x_coordinate_to_xz(x: int) -> ProjectivePointXZ:
    return ProjectivePointXZ(x, K(1))


def xz_to_x_coordinate(P: ProjectivePointXZ) -> int:
    return P.X // P.Z


def xADD(P: ProjectivePointXZ, Q: ProjectivePointXZ, Diff: ProjectivePointXZ) -> ProjectivePointXZ:
    # Montgomery curves and their arithmetic, Algorithm 1
    V0 = P.X + P.Z
    V1 = Q.X - Q.Z
    V1 = V1 * V0
    V0 = P.X - P.Z
    V2 = Q.X + Q.Z
    V2 = V2 * V0
    V3 = V1 + V2
    V3 = V3 ^ 2
    V4 = V1 - V2
    V4 = V4 ^ 2
    X_plus = Diff.Z * V3
    Z_plus = Diff.X * V4

    return ProjectivePointXZ(X_plus, Z_plus)


def xDBL(P: ProjectivePointXZ) -> ProjectivePointXZ:
    # Montgomery curves and their arithmetic, Algorithm 2
    V1 = P.X + P.Z
    V1 = V1 ^ 2
    V2 = P.X - P.Z
    V2 = V2 ^ 2
    X_2P = V1 * V2
    V1 = V1 - V2
    V3 = ((A + 2) / 4) * V1
    V3 = V3 + V2
    Z_2P = V1 * V3

    return ProjectivePointXZ(X_2P, Z_2P)


def SWAP(b: int, x0: int, x1: int):
    # Montgomery curves and their arithmetic, Algorithm 7

    x0 = int(x0)
    x1 = int(x1)

    mask = (b & 1) * 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    v = mask & (x0 ^^ x1)

    return K(int(x0) ^^ int(v)), K(int(x1) ^^ int(v))


def constant_time_montgomery_ladder(P: ProjectivePointXZ, k: int):
    # Montgomery curves and their arithmetic, Algorithm 8
    group_order2 = (
        Integer(
            7237005577332262213973186563042994240857116359379907606001950938285454250989
        )
        * 2
    )
    k = Integer(k) + group_order2
    k_bits = k.bits()

    x0, x1 = xDBL(P), P
    prev_bit = k_bits[k.nbits() - 1]
    for i in range(k.nbits() - 2, -1, -1):
        cur_bit = k_bits[i]

        swap_bit = prev_bit ^^ cur_bit

        x0X, x1X = SWAP(swap_bit, x0.X, x1.X)
        x0Z, x1Z = SWAP(swap_bit, x0.Z, x1.Z)

        x0 = ProjectivePointXZ(x0X, x0Z)
        x1 = ProjectivePointXZ(x1X, x1Z)

        x0, x1 = xDBL(x0), xADD(x0, x1, P)

        prev_bit = cur_bit
    x0X, x1X = SWAP(k_bits[0], x0.X, x1.X)
    x0Z, x1Z = SWAP(k_bits[0], x0.Z, x1.Z)

    x0 = ProjectivePointXZ(x0X, x0Z)
    x1 = ProjectivePointXZ(x1X, x1Z)
    return x0, x1


def y_recovery(G: AffinePoint, Q: ProjectivePointXZ, G_plus_Q: ProjectivePointXZ) -> ProjectivePoint:
    # Montgomery curves and their arithmetic, Algorithm 5
    v1 = G.x * Q.Z
    v2 = Q.X + v1
    v3 = Q.X - v1
    v3 = v3 ^ 2
    v3 = v3 * G_plus_Q.X
    v1 = 2 * A * Q.Z
    v2 = v2 + v1
    v4 = G.x * Q.X
    v4 = v4 + Q.Z
    v2 = v2 * v4
    v1 = v1 * Q.Z
    v2 = v2 - v1
    v2 = v2 * G_plus_Q.Z
    Y_prime = v2 - v3
    v1 = 2 * B * G.y
    v1 = v1 * Q.Z
    v1 = v1 * G_plus_Q.Z
    X_prime = v1 * Q.X
    Z_prime = v1 * Q.Z

    return ProjectivePoint(X_prime, Y_prime, Z_prime)


def montgomery_ladder(P: ProjectivePointXZ, k: int) -> tuple[ProjectivePointXZ, ProjectivePointXZ]:
    # Montgomery curves and their arithmetic, Algorithm 4
    x0, x1 = P, xDBL(P)
    k_bits = Integer(k).bits()
    for i in range(len(k_bits) - 2, -1, -1):
        if k_bits[i] == 0:
            x0, x1 = xDBL(x0), xADD(x0, x1, P)
        else:
            x0, x1 = xADD(x0, x1, P), xDBL(x1)
    return x0, x1
