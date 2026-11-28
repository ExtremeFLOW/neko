"""Compare Euler IDP integration summaries from two execution backends."""

import argparse
import json
import math
from pathlib import Path


FLOAT_TAGS = (
    "EULER_IDP_ERROR",
    "EULER_IDP_DRIFT",
    "EULER_IDP_LIMITER",
    "EULER_IDP_LIMITER_STATS",
    "EULER_IDP_ENTROPY_BLEND",
    "EULER_IDP_STATE",
    "EULER_IDP_BOUNDS",
    "EULER_IDP_CONSERVATION",
    "EULER_IDP_CORRECTION",
    "EULER_IDP_BOUNDARY",
)
INTEGER_TAGS = ("EULER_IDP_LIMITER_COUNTS",)
LIMITER_FLOAT_TAGS = ("EULER_IDP_LIMITER", "EULER_IDP_LIMITER_STATS")


def parse_log(path):
    """Return the final tagged summary in one solver log."""
    tags = {}
    expected = ("EULER_IDP_RESULT", *FLOAT_TAGS, *INTEGER_TAGS)
    for line in path.read_text(errors="replace").splitlines():
        fields = line.strip().split()
        if fields and fields[0] in expected:
            tags[fields[0]] = fields[1:]
    if not all(tag in tags for tag in expected):
        return None
    return tags


def collect_logs(directory):
    """Collect every completed Euler IDP case by log basename."""
    summaries = {}
    for path in sorted(directory.glob("*.log")):
        summary = parse_log(path)
        if summary is not None:
            summaries[path.name] = summary
    return summaries


def tolerances(precision, tag):
    """Return backend-parity tolerances for one summary category."""
    if precision == "sp":
        if tag in LIMITER_FLOAT_TAGS:
            return 3.0e-4, 3.0e-5
        return 1.0e-4, 3.0e-5
    if tag in LIMITER_FLOAT_TAGS:
        return 3.0e-7, 3.0e-9
    return 5.0e-10, 2.0e-11


def integer_tolerance(precision, tag):
    """Return an absolute tolerance for threshold-based integer statistics."""
    if tag == "EULER_IDP_LIMITER_COUNTS":
        return 4 if precision == "sp" else 2
    return 0


def compare(cpu_logs, device_logs, precision):
    """Compare all completed cases common to both log directories."""
    cpu = collect_logs(cpu_logs)
    device = collect_logs(device_logs)
    common = sorted(cpu.keys() & device.keys())
    if not common:
        raise AssertionError("No completed Euler IDP logs were found to compare")
    if cpu.keys() != device.keys():
        missing_from_cpu = sorted(device.keys() - cpu.keys())
        missing_from_device = sorted(cpu.keys() - device.keys())
        raise AssertionError(
            "Completed case sets differ: "
            f"missing from CPU={missing_from_cpu}, "
            f"missing from device={missing_from_device}"
        )

    report = {
        "precision": precision,
        "case_count": len(common),
        "cases": common,
        "metrics": {},
    }
    failures = []
    for name in common:
        cpu_summary = cpu[name]
        device_summary = device[name]
        if cpu_summary["EULER_IDP_RESULT"] != device_summary["EULER_IDP_RESULT"]:
            failures.append(f"{name}: result metadata differs")

        for tag in INTEGER_TAGS:
            cpu_values = [int(value) for value in cpu_summary[tag]]
            device_values = [int(value) for value in device_summary[tag]]
            if len(cpu_values) != len(device_values):
                failures.append(f"{name}: {tag} has different lengths")
                continue
            absolute_tolerance = integer_tolerance(precision, tag)
            metric = report["metrics"].setdefault(
                tag,
                {
                    "absolute_tolerance": absolute_tolerance,
                    "max_absolute_difference": 0,
                },
            )
            for index, (cpu_value, device_value) in enumerate(
                zip(cpu_values, device_values)
            ):
                difference = abs(cpu_value - device_value)
                metric["max_absolute_difference"] = max(
                    metric["max_absolute_difference"], difference
                )
                if difference > absolute_tolerance:
                    failures.append(
                        f"{name}: {tag}[{index}] differs: "
                        f"CPU={cpu_value}, device={device_value}"
                    )

        for tag in FLOAT_TAGS:
            cpu_values = [float(value) for value in cpu_summary[tag]]
            device_values = [float(value) for value in device_summary[tag]]
            if len(cpu_values) != len(device_values):
                failures.append(f"{name}: {tag} has different lengths")
                continue
            relative_tolerance, absolute_tolerance = tolerances(precision, tag)
            metric = report["metrics"].setdefault(
                tag,
                {
                    "absolute_tolerance": absolute_tolerance,
                    "relative_tolerance": relative_tolerance,
                    "max_absolute_difference": 0.0,
                    "max_relative_difference": 0.0,
                },
            )
            for index, (cpu_value, device_value) in enumerate(
                zip(cpu_values, device_values)
            ):
                difference = abs(cpu_value - device_value)
                scale = max(abs(cpu_value), abs(device_value))
                relative_difference = difference / scale if scale else 0.0
                metric["max_absolute_difference"] = max(
                    metric["max_absolute_difference"], difference
                )
                metric["max_relative_difference"] = max(
                    metric["max_relative_difference"], relative_difference
                )
                if not math.isclose(
                    cpu_value,
                    device_value,
                    rel_tol=relative_tolerance,
                    abs_tol=absolute_tolerance,
                ):
                    failures.append(
                        f"{name}: {tag}[{index}] differs: "
                        f"CPU={cpu_value:.17e}, device={device_value:.17e}"
                    )

    if failures:
        raise AssertionError("\n".join(failures))
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu-logs", type=Path, required=True)
    parser.add_argument("--device-logs", type=Path, required=True)
    parser.add_argument("--precision", choices=("sp", "dp"), default="dp")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    report = compare(
        arguments.cpu_logs,
        arguments.device_logs,
        arguments.precision,
    )
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output is not None:
        arguments.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
