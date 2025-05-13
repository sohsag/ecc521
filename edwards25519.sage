from utils import (
    Edwards_Neutral_Element_Affine,
    Edwards_Neutral_Element_Projective,
    AffinePoint,
    ProjectivePoint,
    affine_to_proj,
    proj_to_affine
)
import os
import hashlib
p = 0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffed
K = GF(p)
a = K(0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffec)
d = K(0x52036cee2b6ffe738cc740797779e89800700a4d4141d8ab75eb4dca135978a3)
E = EllipticCurve(K, (K(-1/48) * (a^2 + 14*a*d + d^2),K(1/864) * (a + d) * (-a^2 + 34*a*d - d^2)))
E.set_order(0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed * 0x08)



def to_weierstrass(P: AffinePoint) -> AffinePoint:
	return AffinePoint((5*a + a*P.y - 5*d*P.y - d)/(12 - 12*P.y), (a + a*P.y - d*P.y -d)/(4*P.x - 4*P.x*P.y))

def to_twistededwards(P: AffinePoint) -> AffinePoint:
    u, v = P
    y = (5*a - 12*u - d)/(-12*u - a + 5*d)
    x = (a + a*y - d*y -d)/(4*v - 4*v*y)
    return AffinePoint(x, y)



def edwards_addition(P1: ProjectivePoint, P2: ProjectivePoint) -> ProjectivePoint:
    # Algorithm 3
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
    Y3 = A * G * (D - a * C) # a here
    Z3 = G * F

    return ProjectivePoint(X3, Y3, Z3)

Edwards_Neutral_Element_Projective = ProjectivePoint(K(0), K(1), K(1))

def edwards_scalar_multiplication(P: AffinePoint, scalar: int) -> AffinePoint:
    P_proj = affine_to_proj(P, K)
    Q = P_proj
    R = Edwards_Neutral_Element_Projective
    while scalar > 0:
        if scalar % 2 == 1:
            R = edwards_addition(R, Q)
        Q = edwards_addition(Q, Q)
        scalar = scalar // 2

    return proj_to_affine(R) # er det ikke akavet at være i projective? måske vi burde lave pointen om til projective i funktionen og returnere 

Gx, Gy = K(0x216936D3CD6E53FEC0A4E231FDD6DC5C692CC7609525A7B2C9562D608F25D51A), K(0x6666666666666666666666666666666666666666666666666666666666666658)
B = AffinePoint(Gx, Gy)
l = 2 ^ 252 + 27742317777372353535851937790883648493  # B's group order

def generate_keypair():
    # Implemented accordingly to "High-speed high-security signatures"
    k = os.urandom(32)
    b_least_significant_bits = hashlib.sha512(k).digest()[:32]
    s = int.from_bytes(b_least_significant_bits, 'little')

    A = edwards_scalar_multiplication(B, s)
    return A, k

def encode_point(P: AffinePoint):
    #Implemented according to https://datatracker.ietf.org/doc/html/rfc8032#section-5.1.2
    x, y = P

    y_bytes = list(y.to_bytes('little'))
    if x.lift() & 1: 
        y_bytes[31] |= 0x80
    return bytes(y_bytes)

def decode_point(encoded_bytes):
    # Implemented according to https://datatracker.ietf.org/doc/html/rfc8032#section-5.1.3
    encoded = bytearray(encoded_bytes)
    x_0 = (encoded[31] >> 7) & 1
    encoded[31] &= 0x7F
    y_int = int.from_bytes(encoded, 'little')

    p = 2^255 - 19
    if y_int >= p:
        raise ValueError("Decoding fails: y >= p")

    y = K(y_int)  
    
    u = (y ^ 2 - 1)
    d = -121665 * pow(121666, -1, p) 
    v = (d * y ^ 2 + 1)
    
    pow_result = pow(u * v ^ 7, (p - 5) // 8, p)
    x = (u * v ^ 3 * pow_result)
    
    v_x2 = (v * x ^ 2)
    
    if v_x2 == u:
        pass
    elif v_x2 == (-u):
        x = (x * pow(2, (p - 1) // 4, p))
    else:
        raise ValueError("Decoding fails: no square root exists")
    
    if x == 0 and x_0 == 1:
        raise ValueError("Decoding fails: x = 0 and x_0 = 1")
    
    if x.lift() % 2 != x_0:
        x = p - x
    
    return AffinePoint(x, y)


def sign_message(message, k):
    # Implemented accordingly to "High-speed high-security signatures"
    secret_key_hash = hashlib.sha512(k).digest()
    b_least_significant_bits = secret_key_hash[:32]
    
    a = int.from_bytes(b_least_significant_bits, 'little') % l
    
    A = edwards_scalar_multiplication(B, a)
    A_encoded = encode_point(A)
    message_prefix = secret_key_hash[32:]
    r = hashlib.sha512(message_prefix + message).digest()
    r_int = int.from_bytes(r, 'little') % l

    R = edwards_scalar_multiplication(B, r_int)
    R_encoded = encode_point(R)

    S = (r_int + (int.from_bytes(hashlib.sha512(R_encoded + A_encoded + message).digest(),'little')*a)) % l
    return int.from_bytes(R_encoded, 'little'), S

def verify_message(message, A, R_encoded_int, S):
    # Implemented accordingly to "High-speed high-security signatures"
    R_encoded = int.to_bytes(R_encoded_int, 32, 'little')
    R = decode_point(R_encoded)
    lhs = edwards_scalar_multiplication(edwards_scalar_multiplication(B,S), 8)
    HAM = int.from_bytes(hashlib.sha512(R_encoded + encode_point(A) + message).digest(), 'little') % l
    rhs_1 = edwards_scalar_multiplication(
        edwards_scalar_multiplication(
            A, HAM
        ), 8)
    rhs_2 = edwards_scalar_multiplication(R,8) 
    rhs_1 = affine_to_proj(rhs_1, K)
    rhs_2 = affine_to_proj(rhs_2, K)
    rhs = edwards_addition(rhs_1, rhs_2)
    rhs = proj_to_affine(rhs)
    print(rhs)
    print(lhs)
    if rhs != lhs:
        raise ValueError("Signature not valid")
    
    return "Signature accepted"
