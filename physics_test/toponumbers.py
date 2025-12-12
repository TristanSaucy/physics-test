from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CandidateSet:
    name: str
    values: tuple[float, ...]
    note: str = ""


def candidate_sets() -> list[CandidateSet]:
    """
    Curated candidate C sets (dimensionless integers / simple rationals).

    These are *not* claims of physical meaning — just convenient families to scan.
    """
    def _octaves(base: int, k_min: int, k_max: int) -> tuple[float, ...]:
        # integer-only (k>=0) so results stay integral "topological numbers"
        if k_min < 0:
            k_min = 0
        return tuple(float(base * (2**k)) for k in range(k_min, k_max + 1))

    def _union(*seqs: tuple[float, ...]) -> tuple[float, ...]:
        s: set[float] = set()
        for seq in seqs:
            s.update(seq)
        return tuple(sorted(s))

    # Octave ranges: tuned so we can reach ~1e11..1e12 magnitudes (useful for hitting m≈-128
    # with inverse-gravity-like targets), without being absurdly huge.
    oct_k_min, oct_k_max = 0, 40

    rotation_degrees = (360, 720, 180, 90, 45, 30, 24, 12, 6)
    harmonic_432 = (432, 216, 144, 108, 72, 54, 36, 27, 18, 9)

    rotation_octaves = _union(*(_octaves(int(v), oct_k_min, oct_k_max) for v in rotation_degrees))
    harmonic_432_octaves = _union(*(_octaves(int(v), oct_k_min, oct_k_max) for v in harmonic_432))
    octaves_360 = _octaves(360, oct_k_min, oct_k_max)

    base_sets = [
        CandidateSet(
            "rotation-degrees",
            rotation_degrees,
            "Degree-based rotational symmetry / circle partitions",
        ),
        CandidateSet(
            "rotation-octaves",
            rotation_octaves,
            "Rotation-degree family scaled by powers of two: v*2^k (octaves)",
        ),
        CandidateSet(
            "harmonic-432-family",
            harmonic_432,
            "Commonly cited harmonic integer family (and divisors)",
        ),
        CandidateSet(
            "harmonic-432-octaves",
            harmonic_432_octaves,
            "432 family scaled by powers of two: v*2^k (octaves)",
        ),
        CandidateSet(
            "360-octaves",
            octaves_360,
            "360 scaled by powers of two: 360*2^k",
        ),
        CandidateSet(
            "powers-of-two",
            (512, 256, 128, 64, 32, 16, 8, 4, 2, 1),
            "Powers of two (useful sanity family)",
        ),
        CandidateSet(
            "fibonacci-ish",
            (987, 610, 377, 233, 144, 89, 55, 34, 21, 13, 8, 5, 3, 2, 1),
            "Fibonacci numbers (ties to phi relationships)",
        ),
        CandidateSet(
            "misc-round",
            (1000, 600, 500, 400, 300, 240, 200, 120),
            "Other round-number baselines",
        ),
    ]

    # Add a combined union set for broad octave exploration.
    union = CandidateSet(
        "octave-union",
        tuple(sorted(set(rotation_octaves) | set(harmonic_432_octaves) | set(octaves_360))),
        "Union of octave-scaled families (rotation-octaves + harmonic-432-octaves + 360-octaves)",
    )

    return base_sets + [union]


def get_candidate_set(name: str) -> CandidateSet:
    for s in candidate_sets():
        if s.name == name:
            return s
    raise KeyError(f"Unknown candidate set: {name}")
