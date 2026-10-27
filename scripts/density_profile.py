"""Plot the positive-spin density profile from a saved configuration or density CSV."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Final lattice CSV or averaged density CSV")
    parser.add_argument("--save", type=Path, help="Optional image output path")
    arguments = parser.parse_args()

    if arguments.input.name.endswith("_density.csv"):
        data = pd.read_csv(arguments.input)
        profile = data["mean_plus_density"].to_numpy(dtype=float)
    else:
        lattice = np.loadtxt(arguments.input, delimiter=",")
        if lattice.ndim != 2:
            parser.error("Expected a two-dimensional lattice CSV")
        profile = (lattice == 1).mean(axis=1)

    figure, axis = plt.subplots(figsize=(7, 4))
    axis.plot(np.arange(len(profile)), profile, marker="o")
    axis.set(xlabel="Row y", ylabel="Fraction of +1 spins", ylim=(-0.05, 1.05))
    figure.tight_layout()
    if arguments.save:
        arguments.save.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(arguments.save, dpi=170)
    else:
        plt.show()
    plt.close(figure)


if __name__ == "__main__":
    main()
