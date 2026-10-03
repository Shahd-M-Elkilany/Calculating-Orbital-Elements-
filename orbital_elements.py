"""Orbit determination from three position vectors (Gibbs method).

Given three geocentric position vectors of a satellite taken at successive
times, the Gibbs method recovers the velocity at the middle observation.
The position and velocity are then converted to the six classical orbital
elements.

Units: kilometres, seconds, degrees. Frame: Earth-centred inertial (ECI).
"""
from __future__ import annotations

import argparse

import numpy as np

MU_EARTH = 398600.0  # Earth's gravitational parameter, km^3/s^2


def gibbs_velocity(r1, r2, r3, mu=MU_EARTH, coplanar_tol_deg=3.0):
    """Return the velocity vector (km/s) at r2 using the Gibbs method.

    The three vectors must be coplanar. A ValueError is raised when r1 is
    more than `coplanar_tol_deg` degrees out of the plane of r2 and r3.
    """
    r1, r2, r3 = (np.asarray(r, dtype=float) for r in (r1, r2, r3))
    r1_mag, r2_mag, r3_mag = (np.linalg.norm(r) for r in (r1, r2, r3))

    c12 = np.cross(r1, r2)
    c23 = np.cross(r2, r3)
    c31 = np.cross(r3, r1)

    # Angle between r1 and the normal to the r2-r3 plane; 90 deg when coplanar.
    out_of_plane = np.degrees(
        np.arcsin(abs(np.dot(r1, c23)) / (r1_mag * np.linalg.norm(c23)))
    )
    if out_of_plane > coplanar_tol_deg:
        raise ValueError(
            f"Position vectors are not coplanar ({out_of_plane:.2f} deg out of plane)."
        )

    n = r1_mag * c23 + r2_mag * c31 + r3_mag * c12
    d = c12 + c23 + c31
    s = r1 * (r2_mag - r3_mag) + r2 * (r3_mag - r1_mag) + r3 * (r1_mag - r2_mag)

    n_mag = np.linalg.norm(n)
    d_mag = np.linalg.norm(d)
    return np.sqrt(mu / (n_mag * d_mag)) * (np.cross(d, r2) / r2_mag + s)


def classical_orbital_elements(r, v, mu=MU_EARTH):
    """Convert a state vector to classical orbital elements.

    Returns (a, e, i, raan, argp, nu): semi-major axis in km, eccentricity,
    then inclination, right ascension of the ascending node, argument of
    periapsis and true anomaly in degrees.

    Circular and equatorial orbits leave some angles undefined; those are
    returned as 0.
    """
    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)
    r_mag = np.linalg.norm(r)
    v_mag = np.linalg.norm(v)
    tiny = 1e-10

    h = np.cross(r, v)
    h_mag = np.linalg.norm(h)

    e_vec = ((v_mag**2 - mu / r_mag) * r - np.dot(r, v) * v) / mu
    e = np.linalg.norm(e_vec)

    i = np.arccos(np.clip(h[2] / h_mag, -1.0, 1.0))

    node = np.cross([0.0, 0.0, 1.0], h)
    node_mag = np.linalg.norm(node)

    if node_mag > tiny:
        raan = np.arctan2(node[1], node[0]) % (2 * np.pi)
    else:
        raan = 0.0

    if node_mag > tiny and e > tiny:
        argp = np.arccos(np.clip(np.dot(node, e_vec) / (node_mag * e), -1.0, 1.0))
        if e_vec[2] < 0:
            argp = 2 * np.pi - argp
    else:
        argp = 0.0

    if e > tiny:
        nu = np.arccos(np.clip(np.dot(e_vec, r) / (e * r_mag), -1.0, 1.0))
        if np.dot(r, v) < 0:
            nu = 2 * np.pi - nu
    else:
        nu = 0.0

    energy = v_mag**2 / 2 - mu / r_mag
    a = -mu / (2 * energy)

    return a, e, np.degrees(i), np.degrees(raan), np.degrees(argp), np.degrees(nu)


def _vector(text):
    parts = [float(x) for x in text.split(",")]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("expected three comma-separated values")
    return np.array(parts)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Orbital elements from three position vectors (Gibbs method)."
    )
    parser.add_argument("r1", type=_vector, help="first position, km, as x,y,z")
    parser.add_argument("r2", type=_vector, help="second position, km, as x,y,z")
    parser.add_argument("r3", type=_vector, help="third position, km, as x,y,z")
    parser.add_argument("--mu", type=float, default=MU_EARTH,
                        help="gravitational parameter, km^3/s^2 (default: Earth)")
    args = parser.parse_args(argv)

    v2 = gibbs_velocity(args.r1, args.r2, args.r3, args.mu)
    a, e, i, raan, argp, nu = classical_orbital_elements(args.r2, v2, args.mu)

    print(f"Velocity at r2         : [{v2[0]:.4f}, {v2[1]:.4f}, {v2[2]:.4f}] km/s")
    print(f"Semi-major axis (a)    : {a:.1f} km")
    print(f"Eccentricity (e)       : {e:.4f}")
    print(f"Inclination (i)        : {i:.2f} deg")
    print(f"RAAN (Omega)           : {raan:.2f} deg")
    print(f"Argument of periapsis  : {argp:.2f} deg")
    print(f"True anomaly (nu)      : {nu:.2f} deg")


if __name__ == "__main__":
    main()
