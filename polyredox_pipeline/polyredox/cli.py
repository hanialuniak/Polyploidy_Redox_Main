from __future__ import annotations

import argparse

from . import pipeline, simulate


def main():
    parser = argparse.ArgumentParser(prog="polyredox")
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="Run the analysis from a YAML config")
    run.add_argument("config")
    sim = commands.add_parser("simulate", help="Write a simulated dataset with known answers")
    sim.add_argument("--outdir", required=True)
    sim.add_argument("--scenario", choices=simulate.SCENARIOS, default="effect")
    sim.add_argument("--donors", type=int, default=24)
    sim.add_argument("--nuclei-per-library", type=int, default=250)
    sim.add_argument("--seed", type=int, default=1)
    sim.add_argument("--per-copy-fragments", type=int, default=4000)
    sim.add_argument("--write-fragments", action="store_true")
    args = parser.parse_args()
    if args.command == "run":
        pipeline.run(pipeline.load_config(args.config))
    else:
        simulate.simulate(args.outdir, args.scenario, args.donors, args.nuclei_per_library, args.seed,
                          per_copy_fragments=args.per_copy_fragments, write_fragments=args.write_fragments)
