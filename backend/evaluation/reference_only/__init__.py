"""Reference data for offline evaluation; never imported by app runtime."""

from .reaction_suite import REACTIONS, SPECIES, validate_suite

__all__ = ["REACTIONS", "SPECIES", "validate_suite"]
