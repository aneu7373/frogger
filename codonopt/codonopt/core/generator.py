import random


def choose_random_codon(candidates, seed=None, rng=None):
    """
    Choose a random codon from candidates.

    Reproducibility rules:
      - If rng is provided (random.Random), use it.
      - Else if seed is provided, use a *local* RNG seeded once here (does NOT reseed global RNG).
      - Else fall back to global random.choice.

    This avoids the previous bug where random.seed(seed) was called repeatedly,
    contaminating global RNG state and producing brittle behavior.
    """
    if rng is not None:
        return rng.choice(candidates)
    if seed is not None:
        local = random.Random(seed)
        return local.choice(candidates)
    return random.choice(candidates)
