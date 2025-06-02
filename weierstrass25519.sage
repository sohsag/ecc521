from utils import (
    AffinePoint,
    JacobianPoint,
)

p = 0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFED
K = GF(p)
a = K(19298681539552699237261830834781317975544997444273427339909597334573241639236)
b = K(55751746669818908907645289078257140818241103727901012315294400837956729358436)
E = EllipticCurve([K(0), K(0), K(0), a, b])
Gx = K(19298681539552699237261830834781317975544997444273427339909597334652188435546)
Gy = K(14781619447589544791020593568409986887264606134616475288964881837755586237401)
G: AffinePoint = AffinePoint(Gx, Gy)


def weierstrass_addition(P1: AffinePoint, P2: AffinePoint) -> AffinePoint:
    # Introduction to Mathematical Cryptography Theorem 6.6
    if P1 == "O":
        return P2

    if P2 == "O":
        return P1

    x1, y1 = P1
    x2, y2 = P2

    if x1 == x2 and y1 == -y2:
        return "O"

    if P1 == P2:
        lam = (3 * x1 * x1 + a) / (2 * y1)
    else:
        lam = (y2 - y1) / (x2 - x1)

    x3 = lam^2 - x1 - x2
    y3 = lam * (x1 - x3) - y1

    return AffinePoint(x3, y3)


def weierstrass_scalar_multiplication(P: AffinePoint, scalar: int) -> AffinePoint:
    # Introduction to Mathematical Cryptography table 6.3
    Q = P
    R = "O"
    while scalar > 0:
        if scalar % 2 == 1:
            R = weierstrass_addition(R, Q)
        Q = weierstrass_addition(Q, Q)
        scalar = scalar // 2
    return R

def weierstrass_addition_fermat(P: AffinePoint, Q: AffinePoint) -> AffinePoint:
    # Introduction to Mathematical Cryptography Theorem 6.6
    if P == "O":
        return Q

    if Q == "O":
        return P

    x1, y1 = P
    x2, y2 = Q

    if x1 == x2 and y1 == -y2:
        return "O"

    if P == Q:
        lam = (3 * x1.lift() * x1.lift() + a.lift()) \
         * pow((2 * y1.lift()), p-2, p) % p
    else:
        lam = (y2.lift() - y1.lift()) \
         * pow((x2.lift()-x1.lift()), p-2, p) % p

    x3 = K(lam)^2 - x1 - x2
    y3 = K(lam) * (x1 - x3) - y1

    return AffinePoint(x3, y3)

def weierstrass_scalar_multiplication_fermat(P: AffinePoint, scalar: int) -> AffinePoint:
    # Introduction to Mathematical Cryptography table 6.3
    Q = P
    R = "O"
    while scalar > 0:
        if scalar % 2 == 1:
            R = weierstrass_addition_fermat(R, Q)
        Q = weierstrass_addition_fermat(Q, Q)
        scalar = scalar // 2
    return R

identity_element = JacobianPoint(0, 1, 0)

def weierstrass_jacob_addition(P: JacobianPoint, Q: JacobianPoint) -> JacobianPoint:
    # https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
    X1, Y1, Z1 = P
    X2, Y2, Z2 = Q

    if Z1 == 0:
        return Q
    
    if Z2 == 0:
        return P

    Z1Z1 = Z1^2
    Z2Z2 = Z2^2
    U1 = X1*Z2Z2
    U2 = X2*Z1Z1
    S1 = Y1*Z2*Z2Z2
    S2 = Y2*Z1*Z1Z1

    if U1 == U2:
        if S1 == S2:
            return weierstrass_jacob_double(P)

    H = U2-U1
    I = (2*H)^2
    J = H*I
    r = 2*(S2-S1)
    V = U1*I
    X3 = r^2-J-2*V
    Y3 = r*(V-X3)-2*S1*J
    Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2)*H

    return JacobianPoint(X3, Y3, Z3)

def weierstrass_jacob_double(P: JacobianPoint):
    # https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-2007-bl
    X, Y, Z = P
    XX = X^2
    YY = Y^2
    YYYY = YY^2
    ZZ = Z^2
    S = 2*((X+YY)^2-XX-YYYY)
    M = 3*XX+a*ZZ^2
    T = M^2-2*S
    X2 = T
    Y2 = M*(S-T)-8*YYYY
    Z2 = (Y+Z)^2-YY-ZZ
    return JacobianPoint(X2, Y2, Z2)

def weierstrass_jacob_scalar_multiplication(P: JacobianPoint, scalar: int) -> JacobianPoint:
    # Introduction to Mathematical Cryptography table 6.3
    Q = P
    R = identity_element
    while scalar > 0:
        if scalar % 2 == 1:
            R = weierstrass_jacob_addition(R, Q)
        Q = weierstrass_jacob_addition(Q, Q)
        scalar = scalar // 2
    return R