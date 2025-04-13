from utils import (
    ProjectivePoint,
    AffinePoint,
    ProjectivePointXZ
)
p = 0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed
K = GF(p)
A = K(0x76d06)
B = K(0x01)
E = EllipticCurve(K, ((3 - A^2)/(3 * B^2), (2 * A^3 - 9 * A)/(27 * B^3)))
# E.set_order(0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed * 0x08)
def to_weierstrass(P: AffinePoint) -> AffinePoint:
	return AffinePoint(P.x/B + A/(3*B), P.y/B)

def to_montgomery(P: AffinePoint) -> AffinePoint:
	return AffinePoint(B * (P.x - A/(3*B)), B*P.y)

def normal_proj_to_xz(P: ProjectivePoint) -> ProjectivePointXZ:
    return ProjectivePointXZ(P.X, P.Z)

def x_in_weierstrass(x: int) -> int:
    return x/B + A/(3*B)

def x_coordinate_to_xz(x: int) -> ProjectivePointXZ:
    return ProjectivePointXZ(x, K(1))

def xz_to_x_coordinate(P: ProjectivePointXZ) -> int:
    return P.X // P.Z

def diff(P: ProjectivePointXZ, Q: ProjectivePointXZ) -> ProjectivePointXZ:
    return ProjectivePointXZ(P.X - Q.X, P.Z - Q.Z) # this is wrong. we do this by P + (-Q)

def xADD(P: ProjectivePointXZ, Q: ProjectivePointXZ, Diff: ProjectivePointXZ):
    # Algorithm 1
    V0 = P.X + P.Z
    V1 = Q.X - Q.Z
    V1 = V1 * V0
    V0 = P.X - P.Z
    V2 = Q.X + Q.Z
    V2 = V2 * V0
    V3 = V1 + V2
    V3 = V3^2
    V4 = V1 - V2
    V4 = V4^2
    X_plus = Diff.Z * V3
    Z_plus = Diff.X * V4

    return ProjectivePointXZ(X_plus, Z_plus)

def xDBL(P: ProjectivePointXZ) -> ProjectivePointXZ:
    # Algorithm 2
    V1 = P.X + P.Z
    V1 = V1^2
    V2 = P.X - P.Z
    V2 = V2^2
    X_2P = V1 * V2
    V1 = V1 - V2
    V3 = ((A + 2)/4) * V1
    V3 = V3 + V2
    Z_2P = V1 * V3
    
    return ProjectivePointXZ(X_2P, Z_2P)

def SWAP_coordinate(b: int, x0: int, x1: int):
    # Algorithm 7
    # we need to swap first the x coordinate and then the z coordinate. also pad them (technically we don't need to but for sake of showing we do anyway
    # You're effectively computing [0]P + [k]P, which just equals [k]P

    x0 = int(x0)
    x1 = int(x1)

    mask = 0
    for _ in range(K.characteristic().nbits()-1):
        mask = (mask << 1) | b

    v = mask & (x0 ^^ x1)
    
    return K(int(x0) ^^ int(v)), K(int(x1) ^^ int(v))

# we need to "pad" the integer k with 0, but how we want to do this is take the order of E, 
# and start backwards where we start from log_2(order.E) down to 0 and we shift the bit with the index e.g.

"""
for i in range(log_2(E.order(), -1, -1):
    bit_i = k >> i + 1 & 1
    swap
    # do operations here
    swap


    # also try to do the if the bits change
    swap_bit = k_bits[i+1] ^^ k_bits[i]

    





def uniform_ladder(P: ProjectivePointXZ, k: int):
    # Algorithm 8
    k_bits = Integer(k).bits()
    k_bits.reverse() 
    
    x0 = xDBL(P)
    x1 = P
    
    
    for i in range(len(k_bits)-2, -1, -1):
        swap_bit = k_bits[i+1] ^^ k_bits[i]
        
        (x0, x1) = SWAP_coordinate(swap_bit, x0, x1)
        
        x0 = xDBL(x0)
        x1 = xADD(x0, x1, P)
    
    (x0, x1) = SWAP_coordinate(k_bits[0], x0, x1)
    
    return x0
"""
def ladder(P: ProjectivePointXZ, k: int):
    x0, x1 = xDBL(P), P

    n_bits = ceil(log(E.order(), 2))
    
    for i in range(n_bits-2, -1, -1):
        swap_bit = (k >> i) & 1
    
        x0X, x1X = SWAP_coordinate(swap_bit, x0.X, x1.X)
        x0Z, x1Z = SWAP_coordinate(swap_bit, x0.Z, x1.Z)

        x0 = ProjectivePointXZ(x0X, x0Z)
        x1 = ProjectivePointXZ(x1X, x1Z)

        x0, x1 = xDBL(x0), xADD(x0, x1, P)
        
        x0X, x1X = SWAP_coordinate(swap_bit, x0.X, x1.X)
        x0Z, x1Z = SWAP_coordinate(swap_bit, x0.Z, x1.Z)
        
        x0 = ProjectivePointXZ(x0X, x0Z)
        x1 = ProjectivePointXZ(x1X, x1Z)    
    
    return x0


