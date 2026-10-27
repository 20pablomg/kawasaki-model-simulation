"""Visualize simulation output without external datasets."""
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--save", type=Path)
    args = p.parse_args()
    df = pd.read_csv(args.summary)
    if args.save:
        args.save.mkdir(parents=True, exist_ok=True)
    for metric, ylabel in [("energy_per_site", "Mean energy per spin"),
                            ("specific_heat", "Heat capacity per spin"),
                            ("domain_order", "Domain order parameter"),
                            ("domain_response_proxy", "Domain fluctuation response proxy")]:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for (m,L), group in df.groupby(["m0_requested", "L"]):
            v = group.groupby("T")[metric].agg(["mean", "sem"]).sort_index()
            ax.errorbar(v.index, v["mean"], yerr=v["sem"].fillna(0), marker="o", capsize=3, label=f"L={L}, m0={m:g}")
        ax.set(xlabel="Temperature (J/kB)", ylabel=ylabel)
        ax.legend(fontsize=8)
        fig.tight_layout()
        if args.save: fig.savefig(args.save / f"{metric}.png", dpi=170)
        else: plt.show()
        plt.close(fig)
if __name__ == "__main__": main()
