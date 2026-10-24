"""Create a minimal 12-point Fibonacci sphere in GALAMOST XML format."""

from argparse import ArgumentParser
from pathlib import Path

from shape_model import Building, FibonacciSphere


def build(output_path):
    """Build the example and return the generated XML path."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sphere = FibonacciSphere()
    sphere.setSphere(r0=0.25, R=1.0, tag="A", num_points=12)

    building = Building(lx=3.0, ly=3.0, lz=3.0)
    building.Pack([sphere], filename=str(output_path))
    return output_path


def main(argv=None):
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("fibonacci_sphere.xml"),
        help="output XML path (default: fibonacci_sphere.xml)",
    )
    args = parser.parse_args(argv)

    output_path = build(args.out)
    print(f"Wrote 12 particles to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
