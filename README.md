# Orbit Determination with the Gibbs Method

Computes a satellite's orbit from three position vectors. The Gibbs method
recovers the velocity at the middle observation, and the resulting state
vector is converted to the six classical orbital elements.

## Usage

```bash
pip install numpy
python orbital_elements.py -- "-294.32,4265.1,5986.7" "-1365.5,3637.6,6346.8" "-2940.3,2473.7,6555.8"
```

The `--` keeps the leading minus signs from being read as options.

```
Velocity at r2         : [-6.2174, -4.0122, 1.5990] km/s
Semi-major axis (a)    : 8001.4 km
Eccentricity (e)       : 0.1001
Inclination (i)        : 60.00 deg
RAAN (Omega)           : 40.00 deg
Argument of periapsis  : 30.07 deg
True anomaly (nu)      : 49.93 deg
```

Or from Python:

```python
from orbital_elements import gibbs_velocity, classical_orbital_elements

v2 = gibbs_velocity(r1, r2, r3)
a, e, i, raan, argp, nu = classical_orbital_elements(r2, v2)
```

## Conventions

| Quantity | Convention |
| --- | --- |
| Inputs | Three position vectors in km, in time order |
| Frame | Earth-centred inertial (ECI) |
| Gravitational parameter | 398600 km³/s² (Earth); override with `--mu` |
| Angles | Degrees |
| Output | a, e, i, Ω, ω, ν at the second observation |

The three vectors must lie in one plane. Input more than 3° out of plane is
rejected. For circular or equatorial orbits the undefined angles are
returned as 0.

## Validation

The example above is Example 5.1 from Curtis, *Orbital Mechanics for
Engineering Students*. The reference answer is a = 8000 km, e = 0.1,
i = 60°, Ω = 40°, ω = 30°, ν = 50°. `test_orbital_elements.py` checks the
code against it:

```bash
pip install pytest
pytest
```

## Method

With position magnitudes r₁, r₂, r₃:

```
N = r₁ (r₂ × r₃) + r₂ (r₃ × r₁) + r₃ (r₁ × r₂)
D = r₁ × r₂ + r₂ × r₃ + r₃ × r₁
S = r₁ (r₂ − r₃) + r₂ (r₃ − r₁) + r₃ (r₁ − r₂)

v₂ = sqrt(μ / (N·D)) · (D × r₂ / r₂ + S)
```

The Gibbs method works best when the observations are spread widely along
the orbit. For closely spaced observations the Herrick-Gibbs method is the
usual choice; it is not implemented here.
