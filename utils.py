from collections import namedtuple

AffinePoint = namedtuple("AffinePoint", ["x", "y"])
ProjectivePoint = namedtuple("ProjectivePoint", ["X", "Y", "Z"])
JacobianPoint = namedtuple("JacobianPoint", ["X", "Y", "Z"])
ProjectivePointXZ = namedtuple("ProjectivePointXZ", ["X", "Z"])

Edwards_Neutral_Element_Affine = AffinePoint(0, 1)
Edwards_Neutral_Element_Projective = ProjectivePoint(0, 1, 1)

Montgomery_Neutral_Element_Projective = ProjectivePoint(0, 1, 0)


def affine_to_proj(P: AffinePoint, K) -> ProjectivePoint:
    return ProjectivePoint(P.x, P.y, K(1))


def proj_to_affine(P: ProjectivePoint) -> AffinePoint:
    return AffinePoint(P.X / P.Z, P.Y / P.Z)

def affine_to_jacob(P: AffinePoint, K) -> JacobianPoint:
    return JacobianPoint(P.x, P.y, K(1))

def jacob_to_affine(P: JacobianPoint) -> AffinePoint:
    return AffinePoint(P.X/P.Z**2, P.Y/P.Z**3)
