"""Run a parameter sweep of the compiled Kawasaki executable and combine its summaries."""

import argparse
import csv
from pathlib import Path
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=Path("build/kawasaki"))
    parser.add_argument("--sizes", nargs="+", type=int, default=[16, 24, 32])
    parser.add_argument("--temperatures", nargs="+", type=float, default=[1.5, 2.0, 2.2, 2.3, 2.4, 2.6, 3.0])
    parser.add_argument("--magnetizations", nargs="+", type=float, default=[0.0, 0.25])
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--sweeps", type=int, default=5000)
    parser.add_argument("--burn", type=int, default=1000)
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument("--periodic-y", choices=[0, 1], type=int, default=0)
    parser.add_argument("--out", type=Path, default=Path("outputs/sweep"))
    args = parser.parse_args()

    if args.replicates <= 0 or args.sweeps <= args.burn or args.burn < 0 or args.stride <= 0:
        parser.error("Require replicates>0, sweeps>burn>=0, stride>0")
    args.out.mkdir(parents=True, exist_ok=True)
    records = []
    for mi, magnetization in enumerate(args.magnetizations):
        for length in args.sizes:
            for ti, temperature in enumerate(args.temperatures):
                for replicate in range(args.replicates):
                    seed = 20241000 + mi * 100000 + length * 1000 + ti * args.replicates + replicate
                    prefix = args.out / f"m{magnetization:+.2f}_L{length}_T{temperature:.3f}_r{replicate}"
                    command = [str(args.exe), "--L", str(length), "--T", str(temperature),
                               "--m0", str(magnetization), "--sweeps", str(args.sweeps),
                               "--burn", str(args.burn), "--stride", str(args.stride),
                               "--seed", str(seed), "--periodic-y", str(args.periodic_y),
                               "--output", str(prefix)]
                    subprocess.run(command, check=True)
                    with Path(f"{prefix}_summary.csv").open(newline="", encoding="utf-8") as stream:
                        records.append(next(csv.DictReader(stream)))
    with (args.out / "summary_all.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    main()
