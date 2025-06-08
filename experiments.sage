from time import time
from utils import *

def func_format(f):
    print("="*f.__name__.__len__())
    print(" ".join(f.__name__.split("_")))
    f()
    print("="*f.__name__.__len__())

def main():
    func_format(edwards521_vs_edwards25519)
    func_format(edwards25519_vs_montgomery25519_vs_weierstrass25519)
    func_format(edwards25519_vs_montgomery25519_vs_weierstrass25519_fermat)
    func_format(weierstrass25519_affine_vs_weierstrass25519_jacobian)
    func_format(constant_time_montgomery_ladder_vs_varying_time)
    

def edwards521_vs_edwards25519():
    from edwards521 import edwards_scalar_multiplication, K
    Gx = K(0x752cb45c48648b189df90cb2296b2878a3bfd9f42fc6c818ec8bf3c9c0c6203913f6ecc5ccc72434b1ae949d568fc99c6059d0fb13364838aa302a940a2f19ba6c)
    Gy = K(0x0c)
    G_in_edwards521 = AffinePoint(Gx, Gy)
    t1 = time()
    for _ in range(1000):
        Q = edwards_scalar_multiplication(G_in_edwards521, 1716199415032652428745475199770348304317358825035826352348615864796385795849413675475876651663657849636693659065234142604319282948702542317993421293670108523 + 1)
    t2 = time()
    print(f"edwards521 took on avg. {(t2-t1)/1000:.6f} seconds")

    from edwards25519 import edwards_scalar_multiplication, K
    Gx, Gy = K(0x216936D3CD6E53FEC0A4E231FDD6DC5C692CC7609525A7B2C9562D608F25D51A), K(0x6666666666666666666666666666666666666666666666666666666666666658)
    G_in_edwards25519 = AffinePoint(Gx, Gy)
    t1 = time()
    for _ in range(1000):
        Q = edwards_scalar_multiplication(G_in_edwards25519, 57896044618658097711785492504343953926856930875039260848015607506283634007912 + 1)
    t2 = time()
    print(f"edwards25519 took on avg. {(t2-t1)/1000:.6f} seconds")
    
def edwards25519_vs_montgomery25519_vs_weierstrass25519():
    from edwards25519 import edwards_scalar_multiplication, K
    Gx, Gy = K(0x216936D3CD6E53FEC0A4E231FDD6DC5C692CC7609525A7B2C9562D608F25D51A), K(0x6666666666666666666666666666666666666666666666666666666666666658 + 1)
    G_in_edwards25519 = AffinePoint(Gx, Gy)
    t1 = time()
    for _ in range(1000):
        Q = edwards_scalar_multiplication(G_in_edwards25519, 7237005577332262213973186563042994240857116359379907606001950938285454250989 + 1)
    t2 = time()
    print(f"edwards25519 took on avg. {(t2-t1)/1000:.6f} seconds")

    from montgomery25519 import montgomery_ladder, K, normal_proj_to_xz, y_recovery
    Gx = K(0x09)
    Gy = K(0x20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9)
    G_in_montgomery = AffinePoint(Gx, Gy)
    t1 = time()
    for _ in range(1000):
        G = normal_proj_to_xz(affine_to_proj(G_in_montgomery, K))
        x0, x1 = montgomery_ladder(G, 7237005577332262213973186563042994240857116359379907606001950938285454250989 + 1)
        proj_G = y_recovery(G_in_montgomery, x0, x1)
        proj_to_affine(proj_G)

    t2 = time()
    print(f"montgomery25519 took on avg. {(t2-t1)/1000:.6f} seconds")
    from weierstrass25519 import weierstrass_scalar_multiplication, K, G
    t1 = time()
    for _ in range(1000):
        Q = weierstrass_scalar_multiplication(G, 7237005577332262213973186563042994240857116359379907606001950938285454250989 + 1)
    t2 = time()
    print(f"weierstrass25519 with affine coordinates took on avg. {(t2-t1)/1000:.6f} seconds")
    
def edwards25519_vs_montgomery25519_vs_weierstrass25519_fermat():
    from edwards25519 import edwards_scalar_multiplication, K
    Gx, Gy = K(0x216936D3CD6E53FEC0A4E231FDD6DC5C692CC7609525A7B2C9562D608F25D51A), K(0x6666666666666666666666666666666666666666666666666666666666666658)
    G_in_edwards25519 = AffinePoint(Gx, Gy)
    t1 = time()
    for _ in range(1000):
        Q = edwards_scalar_multiplication(G_in_edwards25519, 7237005577332262213973186563042994240857116359379907606001950938285454250989 + 1)
    t2 = time()
    print(f"edwards25519 took on avg. {(t2-t1)/1000:.6f} seconds")

    from montgomery25519 import montgomery_ladder, K, normal_proj_to_xz, y_recovery
    Gx = K(0x09)
    Gy = K(0x20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9)
    G_in_montgomery = AffinePoint(Gx, Gy)
    G = normal_proj_to_xz(affine_to_proj(G_in_montgomery, K))
    t1 = time()
    for _ in range(1000):
        G = normal_proj_to_xz(affine_to_proj(G_in_montgomery, K))
        x0, x1 = montgomery_ladder(G, 7237005577332262213973186563042994240857116359379907606001950938285454250989 + 1)
        proj_G = y_recovery(G_in_montgomery, x0, x1)
        proj_to_affine(proj_G)
    t2 = time()
    print(f"montgomery25519 took on avg. {(t2-t1)/1000:.6f} seconds")

    from weierstrass25519 import weierstrass_scalar_multiplication_fermat, K, G
    t1 = time()
    for _ in range(1000):
        Q = weierstrass_scalar_multiplication_fermat(G, 7237005577332262213973186563042994240857116359379907606001950938285454250989 + 1)
    t2 = time()
    print(f"weierstrass25519 with fermat and affine coordinates took on avg. {(t2-t1)/1000:.6f} seconds")
    
def constant_time_montgomery_ladder_vs_varying_time():
    from montgomery25519 import montgomery_ladder, K, normal_proj_to_xz, y_recovery, constant_time_montgomery_ladder
    Gx = K(0x09)
    Gy = K(0x20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9)
    G_in_montgomery = AffinePoint(Gx, Gy)
    G = normal_proj_to_xz(affine_to_proj(G_in_montgomery, K))
    t1 = time()
    for _ in range(1000):
        G = normal_proj_to_xz(affine_to_proj(G_in_montgomery, K))
        x0, x1 = constant_time_montgomery_ladder(G, 7237005577332262213973186563042994240857116359379907606001950938285454250989 + 1)
        proj_G = y_recovery(G_in_montgomery, x0, x1)
        proj_to_affine(proj_G)
    t2 = time()
    print(f"montgomery25519 constant time took on avg. {(t2-t1)/1000:.6f} seconds")

    Gx = K(0x09)
    Gy = K(0x20ae19a1b8a086b4e01edd2c7748d14c923d4d7e6d7c61b229e9c5a27eced3d9)
    G_in_montgomery = AffinePoint(Gx, Gy)
    G = normal_proj_to_xz(affine_to_proj(G_in_montgomery, K))
    t1 = time()
    for _ in range(1000):
        G = normal_proj_to_xz(affine_to_proj(G_in_montgomery, K))
        x0, x1 = montgomery_ladder(G, 7237005577332262213973186563042994240857116359379907606001950938285454250989 + 1)
        proj_G = y_recovery(G_in_montgomery, x0, x1)
        proj_to_affine(proj_G)
    t2 = time()
    print(f"montgomery25519 varying time took on avg. {(t2-t1)/1000:.6f} seconds")

def weierstrass25519_affine_vs_weierstrass25519_jacobian():
    from weierstrass25519 import weierstrass_scalar_multiplication_fermat, weierstrass_jacob_scalar_multiplication, K, G
    t1 = time()
    for _ in range(1000):
        Q = weierstrass_scalar_multiplication_fermat(G, 7237005577332262213973186563042994240857116359379907606001950938285454250989 + 1)
    t2 = time()
    print(f"weierstrass25519 with affine coordinates took on avg. {(t2-t1)/1000:.6f} seconds")

    t1 = time()
    for _ in range(1000):
        G_proj = affine_to_jacob(G, K)
        Q = weierstrass_jacob_scalar_multiplication(G_proj, 7237005577332262213973186563042994240857116359379907606001950938285454250989 + 1)
        _ = jacob_to_affine(Q)
    t2 = time()
    print(f"weierstrass25519 with jacobian coordinates took on avg. {(t2-t1)/1000:.6f} seconds")


main()
