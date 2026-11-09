from __future__ import annotations

import argparse
import csv
import random
import sys
import time
from pathlib import Path
import networkx as nx
HERE = Path(__file__).resolve().parent
ALGORITHM_DIR = HERE / "algorithm"
sys.path.insert(0, str(ALGORITHM_DIR))

import algorithm.baselines as baselines
import algorithm.our_algorithms as algorithms

for module in (algorithms, baselines):
    if Path(module.__file__).resolve().parent != ALGORITHM_DIR.resolve():
        raise ImportError(f"algorithm loaded from unexpected path: {module.__file__}")

PIN = algorithms.PIN
P_RSR = algorithms.P_RSR
P_IR = algorithms.P_IR
P_BIR = algorithms.P_BIR
ADR = baselines.ADR
TRPW = baselines.TRPW
ABM = baselines.ABM

GRAPH_PATH = HERE / "dataset" / "email-Enron" / "email-Enron-p.txt"
CNEXP_PATH = HERE / "dataset" / "email-Enron" / "email-Enron_cnexp.dat"
OUTPUT_DIR = HERE / "output"

ALL_ALGORITHMS = ("pin", "p_rsr", "p_ir", "p_bir", "adr", "trpw", "abm")
DISPLAY_NAMES = {
    "pin": "PIN",
    "p_rsr": "P+RSR",
    "p_ir": "P+IR",
    "p_bir": "P+BIR",
    "adr": "ADR",
    "trpw": "TRPW",
    "abm": "ABM",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run all representative-world algorithms on one uncertain graph."
    )
    parser.add_argument(
        "--graph",
        type=Path,
        default=GRAPH_PATH,
        help=f"uncertain graph file (default: {GRAPH_PATH})",
    )
    parser.add_argument(
        "--cnexp",
        type=Path,
        default=CNEXP_PATH,
        help=f"CN-expectation file (default: {CNEXP_PATH})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help=f"output directory (default: {OUTPUT_DIR})",
    )
    parser.add_argument(
        "--algorithms",
        nargs="+",
        choices=ALL_ALGORITHMS,
        default=list(ALL_ALGORITHMS),
        help="algorithms to run (default: all)",
    )
    parser.add_argument("--seed", type=int, default=42, help="random seed")
    parser.add_argument(
        "--separator", default=",", help="input field separator (default: ',')"
    )

    parser.add_argument(
        "--p-rsr-factor",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "--p-ir-factor",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "--adr-factor",
        type=float,
        default=10.0,
    )
    parser.add_argument(
        "--trpw-factor",
        type=float,
        default=10.0,
    )
    parser.add_argument(
        "--th-q", type=float, default=0.02, 
    )
    parser.add_argument(
        "--th-c", type=float, default=0.9, 
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not args.graph.is_file():
        raise FileNotFoundError(f"graph file does not exist: {args.graph}")
    if args.cnexp is not None and not args.cnexp.is_file():
        raise FileNotFoundError(f"CN-expectation file does not exist: {args.cnexp}")
    if not args.separator:
        raise ValueError("separator must not be empty")
    for name in ("p_rsr_factor", "p_ir_factor", "adr_factor", "trpw_factor"):
        if getattr(args, name) < 0:
            raise ValueError(f"{name.replace('_', '-')} must be non-negative")
    if not 0.0 < args.th_q < 1.0:
        raise ValueError("th-q must be in (0, 1)")
    if not 0.0 < args.th_c < 1.0:
        raise ValueError("th-c must be in (0, 1)")


def _iter_rows(path: Path, separator: str):
    with path.open("r", encoding="utf-8") as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = [field.strip() for field in line.split(separator)]
            if len(fields) != 3:
                raise ValueError(
                    f"{path}:{line_number}: expected 3 fields separated by "
                    f"{separator!r}, got {len(fields)}"
                )
            yield line_number, fields


def load_uncertain_graph(path: Path, separator: str) -> nx.Graph:
    graph = nx.Graph()
    for line_number, (raw_u, raw_v, raw_p) in _iter_rows(path, separator):
        try:
            u, v, probability = int(raw_u), int(raw_v), float(raw_p)
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number}: invalid u,v,p values") from exc
        if not 0.0 <= probability <= 1.0:
            raise ValueError(
                f"{path}:{line_number}: edge probability must be in [0, 1]"
            )
        graph.add_edge(u, v, p=probability)
    if graph.number_of_edges() == 0:
        raise ValueError(f"graph contains no edges: {path}")
    return graph


def load_cnexp(path: Path, separator: str) -> dict[tuple[int, int], float]:
    cnexp: dict[tuple[int, int], float] = {}
    for line_number, (raw_u, raw_v, raw_value) in _iter_rows(path, separator):
        try:
            u, v, value = int(raw_u), int(raw_v), float(raw_value)
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number}: invalid u,v,value fields") from exc
        cnexp[(min(u, v), max(u, v))] = value
    return cnexp


def compute_cnexp(graph: nx.Graph) -> dict[tuple[int, int], float]:
    neighbors = {node: set(graph[node]) for node in graph.nodes()}
    cnexp: dict[tuple[int, int], float] = {}
    for u, v in graph.edges():
        value = sum(
            graph[u][w]["p"] * graph[v][w]["p"]
            for w in neighbors[u] & neighbors[v]
        )
        cnexp[(min(u, v), max(u, v))] = value
    return cnexp


def validate_cnexp(graph: nx.Graph, cnexp: dict[tuple[int, int], float]) -> None:
    missing = [
        (min(u, v), max(u, v))
        for u, v in graph.edges()
        if (min(u, v), max(u, v)) not in cnexp
    ]
    if missing:
        preview = ", ".join(map(str, missing[:5]))
        suffix = " ..." if len(missing) > 5 else ""
        raise ValueError(
            f"CN-expectation file is missing {len(missing)} graph edges: "
            f"{preview}{suffix}"
        )


def compute_mad(
    graph: nx.Graph, world: nx.Graph, cnexp: dict[tuple[int, int], float]
) -> float:
    neighbors = {node: set(world[node]) for node in graph.nodes()}
    total = 0.0
    for u, v in graph.edges():
        key = (min(u, v), max(u, v))
        total += abs(len(neighbors[u] & neighbors[v]) - cnexp[key])
    return total / graph.number_of_edges()


def write_world(world: nx.Graph, path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        for u, v in sorted((min(u, v), max(u, v)) for u, v in world.edges()):
            writer.writerow((u, v))


def write_cnexp(cnexp: dict[tuple[int, int], float], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        for (u, v), value in sorted(cnexp.items()):
            writer.writerow((u, v, f"{value:.17g}"))


def run_algorithm(
    name: str,
    graph: nx.Graph,
    cnexp: dict[tuple[int, int], float],
    args: argparse.Namespace,
) -> tuple[nx.Graph, int | None, str]:
    """Run one algorithm and return (world, iterations, parameter_summary)."""
    n_edges = graph.number_of_edges()
    expected_edges = sum(graph[u][v]["p"] for u, v in graph.edges())

    random.seed(args.seed)

    if name == "pin":
        return PIN(graph), None, f"seed={args.seed}"
    if name == "p_rsr":
        steps = round(n_edges * args.p_rsr_factor)
        world, _ = P_RSR(graph, PIN(graph), k=steps, cnexp=cnexp)
        return world, steps, f"p_rsr_factor={args.p_rsr_factor}"
    if name == "p_ir":
        steps = round(n_edges * args.p_ir_factor)
        world, _ = P_IR(graph, PIN(graph), k=steps, cnexp=cnexp)
        return world, steps, f"p_ir_factor={args.p_ir_factor}"
    if name == "p_bir":
        world, _ = P_BIR(
            graph, PIN(graph), th_q=args.th_q, th_c=args.th_c, cnexp=cnexp
        )
        return world, None, f"th_q={args.th_q};th_c={args.th_c}"
    if name == "adr":
        steps = round(expected_edges * args.adr_factor)
        return ADR(graph, steps=steps), steps, f"adr_factor={args.adr_factor}"
    if name == "trpw":
        steps = round(expected_edges * args.trpw_factor)
        return TRPW(graph, steps=steps), steps, f"trpw_factor={args.trpw_factor}"
    if name == "abm":
        return ABM(graph), None, ""
    raise ValueError(f"unknown algorithm: {name}")


def main() -> int:
    args = parse_args()
    try:
        validate_args(args)
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Using algorithms from: {ALGORITHM_DIR.resolve()}")
    print(f"  our_algorithms: {Path(algorithms.__file__).resolve()}")
    print(f"  baselines:      {Path(baselines.__file__).resolve()}")
    print(f"Loading graph: {args.graph}")
    try:
        graph = load_uncertain_graph(args.graph, args.separator)
        print(
            f"  nodes={graph.number_of_nodes():,}, edges={graph.number_of_edges():,}"
        )

        if args.cnexp is None:
            print("Computing CN expectations ...")
            start = time.perf_counter()
            cnexp = compute_cnexp(graph)
            cnexp_path = output_dir / "computed_cnexp.dat"
            write_cnexp(cnexp, cnexp_path)
            print(f"  done in {time.perf_counter() - start:.3f}s; saved to {cnexp_path}")
        else:
            print(f"Loading CN expectations: {args.cnexp}")
            cnexp = load_cnexp(args.cnexp, args.separator)
        validate_cnexp(graph, cnexp)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    failed = 0

    selected_names = [DISPLAY_NAMES[name] for name in args.algorithms]
    print(f"Running {len(args.algorithms)} algorithms: {', '.join(selected_names)}")
    for index, name in enumerate(args.algorithms, start=1):
        display_name = DISPLAY_NAMES[name]
        print(f"[{index}/{len(args.algorithms)}] {display_name} ...", flush=True)
        start = time.perf_counter()
        try:
            world, _, _ = run_algorithm(name, graph, cnexp, args)
            runtime = time.perf_counter() - start
            mad = compute_mad(graph, world, cnexp)
            output_path = output_dir / f"email_{name}_output.dat"
            write_world(world, output_path)
            print(
                f"  done: edges={world.number_of_edges():,}, mad={mad:.6f}, "
                f"runtime={runtime:.3f}s"
            )
        except Exception as exc:  
            runtime = time.perf_counter() - start
            failed += 1
            print(f"  failed: {type(exc).__name__}: {exc}", file=sys.stderr)

    if failed:
        print(f"Completed with {failed} failed algorithm(s).", file=sys.stderr)
        return 1
    print("All algorithms completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
