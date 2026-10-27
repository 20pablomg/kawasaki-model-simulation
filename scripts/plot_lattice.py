"""Display a two-dimensional spin configuration saved by the Kawasaki simulator."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path, help="CSV produced as *_final.csv")
    parser.add_argument("--save", type=Path, help="Optional image output path")
    args = parser.parse_args()

    lattice = np.loadtxt(args.snapshot, delimiter=",", dtype=int)
    if lattice.ndim != 2 or lattice.shape[0] != lattice.shape[1]:
        parser.error("Expected a square lattice")
    if not np.isin(lattice, (-1, 1)).all():
        parser.error("Lattice contains values other than -1 and +1")

    fig, ax = plt.subplots(figsize=(5.4, 5))
    image = ax.imshow(lattice, cmap="binary", vmin=-1, vmax=1, interpolation="nearest")
    ax.set(xlabel="x", ylabel="y", title=f"Spin configuration (L={lattice.shape[0]})")
    colorbar = fig.colorbar(image, ax=ax, ticks=[-1, 1], shrink=0.76)
    colorbar.set_label("Spin")
    fig.tight_layout()
    if args.save:
        args.save.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.save, dpi=180)
    else:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
