"""Check against Example 5.1 in Curtis, Orbital Mechanics for Engineering Students."""
import numpy as np

from orbital_elements import classical_orbital_elements, gibbs_velocity

R1 = [-294.32, 4265.1, 5986.7]
R2 = [-1365.5, 3637.6, 6346.8]
R3 = [-2940.3, 2473.7, 6555.8]


def test_gibbs_velocity_matches_reference():
    v2 = gibbs_velocity(R1, R2, R3)
    np.testing.assert_allclose(v2, [-6.2174, -4.0122, 1.5990], atol=2e-3)


def test_elements_match_reference():
    v2 = gibbs_velocity(R1, R2, R3)
    a, e, i, raan, argp, nu = classical_orbital_elements(R2, v2)
    assert abs(a - 8000) < 5
    assert abs(e - 0.1) < 1e-3
    assert abs(i - 60) < 0.05
    assert abs(raan - 40) < 0.05
    assert abs(argp - 30) < 0.5
    assert abs(nu - 50) < 0.5


def test_non_coplanar_vectors_are_rejected():
    try:
        gibbs_velocity([7000, 0, 0], [0, 7000, 0], [0, 0, 7000])
    except ValueError:
        return
    raise AssertionError("expected ValueError")
