"""Download and validate the v0.1 MACE checkpoint."""

from ase.build import molecule

from app.models.mace_adapter import MACECalculatorAdapter


def main() -> None:
    adapter = MACECalculatorAdapter(model="small-0b", device="auto")
    energy = adapter.energy(molecule("H2O"))
    print(f"Loaded {adapter.metadata['model']} on {adapter.metadata['device']}")
    print(f"H2O smoke-test energy: {energy:.6f} eV")


if __name__ == "__main__":
    main()
