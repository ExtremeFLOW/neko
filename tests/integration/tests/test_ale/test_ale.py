"""Integration test for Neko's ALE (moving mesh).

The case is a small cylinder in cross-flow whose surface zone (7) is registered
as a rigid ALE body.  The body motion is the superposition of

  * the built-in harmonic translation/rotation configured in the case file, and
  * an extra contribution injected through the ``ale_rigid_kinematics`` user
    hook in ``test_ale.f90``,

so the test exercises the fact that the two are *added*, not that one replaces
the other.

Because the resulting rigid-body motion has a closed form, the pivot position,
pivot velocity and rotation angle are asserted against the analytic solution
rather than against recorded numbers.  Only the fluid forces/torque are
regression values.

The simulation is run twice:

    part 1:  0 -> END_TIME                 (reference, writes checkpoints)
    part 2:  CHECKPOINT_TIME -> END_TIME   (restarted from part 1's checkpoint)

and every metric in the overlapping window is required to agree.

Log formats parsed here:

    <tstep> <time> <total> <pressure> <viscous> , <label>
    Total_Pivot_pos  <tstep>  <time>  <body>  <x>  <y>  <z>
    Total_Pivot_vel  <tstep>  <time>  <body>  <x>  <y>  <z>
    Total_Rot_deg    <tstep>  <time>  <body>  <roll>  <pitch>  <yaw>
"""

import json
import logging
import os
import shutil
import subprocess
from collections import defaultdict, namedtuple
from math import cos, pi, sin
from os.path import abspath, dirname, isfile, join
from pathlib import Path

import json5
import pytest
from numpy.testing import assert_allclose

import conftest
from conftest import logger
from testlib import configure_nprocs, get_makeneko, run_neko

HERE = abspath(dirname(__file__))
USER_FILE = "test_ale.f90"
CASE_FILE = "test_ale.case"
MESH_NAME = "small_test_cyl.nmsh"

# Name of the ALE body in the case file.
BODY_NAME = "cylinder"

#: Base-shape (phi) field written by part 1 when output_base_shape is on. Part 2
#: imports this.
PHI_FILE = f"phi_{BODY_NAME}0.f00000"


def use_phi_from_first_run(cfg):
    """Point the restart run at the phi field the first run wrote.

    Requires output_base_shape on the first run (it is, in the case file) so
    that phi_<body>0.f00000 exists. Sets import_base_shape on the ALE solver
    and names the file on every registered body, as Neko requires both.
    """
    assert isfile(PHI_FILE), (
        f"{PHI_FILE} was not written by the first run; cannot import it on "
        f"restart. Check that ale.solver.output_base_shape is true."
    )
    cfg["fluid"]["ale"]["solver"]["import_base_shape"] = True
    for body in cfg["fluid"]["ale"]["bodies"]:
        body["base_shape_import_file"] = PHI_FILE

# ---------------------------------------------------------------------------
# Case timeline.
# ---------------------------------------------------------------------------
DT = 1.0e-3
END_TIME = 0.010
CHECKPOINT_TIME = 0.005
#: Step offset between the restarted run and the reference run.
RESTART_OFFSET = round(CHECKPOINT_TIME / DT)
#: Minimum number of time steps the two runs must have in common.
MIN_OVERLAP = 3
#: Instant at which the reference forces below were recorded. Changing this
#: invalidates REFERENCE_FORCES.
COMPARE_TIME = 0.009
COMPARE_STEP_P1 = round(COMPARE_TIME / DT)
COMPARE_STEP_P2 = COMPARE_STEP_P1 - RESTART_OFFSET

# ---------------------------------------------------------------------------
# Analytic rigid-body kinematics.
#
# These must mirror the case file *and* the user file. The case file supplies
# ``oscillation`` (amplitude, frequency) and ``rotation`` (amplitude_deg,
# frequency); the user hook adds the second entry of each pair below.
#
#   v(t) = A 2 pi f cos(2 pi f t)   ->   x(t) = x0 + A sin(2 pi f t)
# ---------------------------------------------------------------------------
PIVOT_0 = (3.0, 2.0, 0.0)
#: (amplitude, frequency) per direction: case-file value + user-hook value.
TRANSLATION = {
    "x": (0.10 - 0.15, 0.1),
    "y": (0.30 - 0.35, 0.2),
}
#: (amplitude in degrees, frequency) contributions to the rotation about z.
ROTATION_Z_DEG = [(2.0 + 3.0, 0.1), (2.0, 0.3)]


def analytic_pivot_position(t):
    """Pivot (x, y) at time `t`, in the absence of time-integration error."""
    return tuple(
        x0 + a * sin(2 * pi * f * t)
        for x0, (a, f) in zip(PIVOT_0, TRANSLATION.values())
    )


def analytic_pivot_velocity(t):
    """Pivot velocity (vx, vy) at time `t`."""
    return tuple(
        a * 2 * pi * f * cos(2 * pi * f * t) for a, f in TRANSLATION.values()
    )


def analytic_rotation_z_deg(t):
    """Rotation angle about z at time `t`, in degrees."""
    return sum(a * sin(2 * pi * f * t) for a, f in ROTATION_Z_DEG)


# ---------------------------------------------------------------------------
# Regression values for the fluid force/torque at COMPARE_TIME.
# ---------------------------------------------------------------------------
REFERENCE_FORCES = {
    # Recorded with long_print = true.
    "dp": {
        "forcex": +1.444186326000e+01,
        "forcey": +7.769968778000e+00,
        "torquez": +2.217687332000e+01,
    },
    # Recorded from the part 1 reference run with long_print = true.
    "sp": {
        "forcex": +1.444019318000e+01,
        "forcey": +7.767942905000e+00,
        "torquez": +2.217312050000e+01,
    },
}

#: Recorded kinematics at COMPARE_TIME.
REFERENCE_KINEMATICS = {
    "dp": {
        "px": +2.9997172582e+00,
        "py": +1.9994345253e+00,
        "vx": -3.1415424236e-02,
        "vy": -6.2827834701e-02,
        "rot_z": +6.2201764966e-02,
    },
    "sp": {
        "px": +2.9997167587e+00,
        "py": +1.9994345903e+00,
        "vx": -3.1415432692e-02,
        "vy": -6.2827795744e-02,
        "rot_z": +6.2257144600e-02,
    },
}

# Tolerances per real precision.
#   kin_*      : analytic pivot position/velocity
#   angle_atol : analytic rotation angle, in degrees
#   force_rtol : regression forces
#   restart_*  : part 1 vs part 2 over the overlap window.
TOLERANCES = {
    "dp": {
        "kin_rtol": 1.0e-8,
        "kin_atol": 1.0e-12,
        "angle_atol": 1.0e-6,
        "force_rtol": 1.0e-6,
        "kin_recorded_rtol": 1.0e-9,
        "restart_rtol": 1.0e-9,
        "restart_atol": 1.0e-12,
    },
    "sp": {
        "kin_rtol": 1.0e-5,
        "kin_atol": 1.0e-7,
        "angle_atol": 1.0e-3,
        "force_rtol": 3.0e-4,
        "kin_recorded_rtol": 1.0e-5,
        "restart_rtol": 1.0e-4,
        "restart_atol": 1.0e-7,
    },
}

#: Timeline for the restart-exactness test. It compares the two runs against
#: each other.
FAST_CHECKPOINT_TIME = 0.004
FAST_END_TIME = 0.007
FAST_RESTART_OFFSET = round(FAST_CHECKPOINT_TIME / DT)

#: Solver tolerance for the same test.
FAST_SOLVER_TOL = {"dp": 1.0e-7, "sp": 1.0e-4}

PRESSURE_FIELD = "Pressure"

#: Tolerance for the restart test.
EXACT_RESTART_RTOL = {"dp": 1.0e-7, "sp": 1.0e-4}

#: Allowed difference in Krylov iteration count between the two runs.
MAX_ITER_DIFF = 3

FORCE_LABELS = ("forcex", "forcey", "forcez", "torquex", "torquey", "torquez")
FORCE_KEYS = ("forcex", "forcey", "torquez")
KINEMATIC_KEYS = ("px", "py", "vx", "vy", "rot_z")

#: One row of the Krylov monitor table.
KspResult = namedtuple("KspResult", "iters res_start res_final")

#: Leading tag
ALE_LOG_TAGS = {
    "Total_Rot_deg": ("rot_x", "rot_y", "rot_z"),
    "Total_Pivot_pos": ("px", "py", "pz"),
    "Total_Pivot_vel": ("vx", "vy", "vz"),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def resolve_executable(name, base=None):
    """Absolute path to an executable or script.

    Relative paths are resolved against `base` (the directory pytest was
    started in) rather than the current one, because the tests chdir into a
    temporary directory before these are needed.
    """
    path = Path(name)
    if base is not None and not path.is_absolute():
        candidate = Path(base) / path
        if candidate.is_file():
            return str(candidate.resolve())
    if path.is_file():
        return str(path.resolve())
    found = shutil.which(name)
    if found is None:
        pytest.fail(
            f"Could not locate '{name}'. Put neko/makeneko on $PATH, or set "
            f"NEKO_EXEC / MAKENEKO_EXEC, or pass an absolute "
            f"--launcher-script."
        )
    return str(Path(found).resolve())


def find_mesh(name):
    """Locate `meshes/<name>`, walking up from this file.

    Independent of the current working directory, so the test keeps working if
    it is moved into a subdirectory of tests/integration.
    """
    directory = HERE
    for _ in range(4):
        candidate = join(directory, "meshes", name)
        if isfile(candidate):
            return candidate
        directory = dirname(directory)
    pytest.fail(f"Could not find meshes/{name} above {HERE}")


def tail(path, n_lines=40):
    """Last `n_lines` of a file, for failure messages."""
    if not isfile(path):
        return f"(no log at {path})"
    with open(path, "r", errors="replace") as f:
        lines = f.readlines()
    return "".join(lines[-n_lines:])


def parse_ale_log(log_path, body=BODY_NAME):
    """Parse a Neko log into ``{tstep: {metric: value}}``.

    Picks up the force/torque lines written by the `force_torque` simulation
    component and the pivot/rotation lines written by `neko_ale%log_pivot` and
    `neko_ale%log_rot_angles`. Only lines belonging to `body` are kept.
    """
    if not isfile(log_path):
        pytest.fail(f"Neko produced no log at {log_path}")

    steps = defaultdict(dict)
    with open(log_path, "r", errors="replace") as f:
        for line in f:
            tokens = line.replace(",", " ").split()
            if len(tokens) < 6:
                continue

            tag = tokens[0]
            if tag in ALE_LOG_TAGS:
                # <tag> <tstep> <time> <body> <x> <y> <z>
                if len(tokens) != 7 or tokens[3] != body:
                    continue
                try:
                    step, time = int(tokens[1]), float(tokens[2])
                    values = [float(v) for v in tokens[4:7]]
                except ValueError:
                    continue
                steps[step].update(zip(ALE_LOG_TAGS[tag], values))
                # The force lines carry `time` at higher precision.
                steps[step].setdefault("t", time)

            elif len(tokens) == 6 and tokens[-1].lower() in FORCE_LABELS:
                # <tstep> <time> <total> <pressure> <viscous> , <label>
                try:
                    step, time = int(tokens[0]), float(tokens[1])
                    total = float(tokens[2])
                except ValueError:
                    continue
                steps[step][tokens[-1].lower()] = total
                steps[step]["t"] = time

    return dict(steps)


def parse_residuals(log_path):
    """Parse the Krylov monitor table into ``{field: {tstep: KspResult}}``.

    Line format, from ``krylov_monitor_print_result``:

        <tstep> | <field name> <iters> <res_start> <res_final>
    """
    residuals = defaultdict(dict)
    with open(log_path, "r", errors="replace") as f:
        for line in f:
            if "|" not in line:
                continue
            head, _, rest = line.partition("|")
            tokens = rest.split()
            if len(tokens) < 4:
                continue
            try:
                step = int(head.strip())
                result = KspResult(
                    iters=int(tokens[-3]),
                    res_start=float(tokens[-2]),
                    res_final=float(tokens[-1]),
                )
            except ValueError:
                continue  # header line, or some other pipe-containing line
            residuals[" ".join(tokens[:-3])][step] = result
    return dict(residuals)


def get_field_residuals(residuals, field, log_path):
    """Residuals for one solver field, failing with the available names."""
    if field not in residuals:
        pytest.fail(
            f"No '{field}' entries in the Krylov monitor table of {log_path}. "
            f"Fields found: {sorted(residuals)}"
        )
    return residuals[field]


def overlap_steps(reference, restarted, minimum=MIN_OVERLAP,
                  offset=RESTART_OFFSET):
    """Steps present in both runs, in the restarted run's numbering."""
    overlap = sorted(s for s in restarted if s + offset in reference)
    assert len(overlap) >= minimum, (
        f"Only {len(overlap)} common time steps between the two runs "
        f"(reference has {sorted(reference)}, restarted has {sorted(restarted)}). "
        f"Expected at least {minimum}."
    )
    return overlap


def require_metrics(parsed, step, keys, log_path, what):
    """Fail with a useful message if `step` is missing any of `keys`."""
    if step not in parsed:
        pytest.fail(
            f"{what}: no metrics logged at time step {step}. Neko probably did "
            f"not get that far.\nLast lines of {log_path}:\n{tail(log_path)}"
        )
    missing = [k for k in keys if k not in parsed[step]]
    if missing:
        pytest.fail(
            f"{what}: time step {step} is missing {missing}. Got "
            f"{sorted(parsed[step])}.\nLast lines of {log_path}:\n"
            f"{tail(log_path)}"
        )


def assert_ale_link_active(log_path, expected="body_attached"):
    """Check that force_torque really linked to the ALE body.

    `setup_ale_link` silently falls back to a fixed centre when the zone is not
    matched to an ALE body, logging "fixed (reverted from ...)".
    """
    with open(log_path, "r", errors="replace") as f:
        centre_lines = [ln.strip() for ln in f if "Center Type:" in ln]

    assert centre_lines, (
        f"No 'Center Type:' line in {log_path}; the force_torque component did "
        f"not initialise."
    )
    for line in centre_lines:
        assert "reverted" not in line, (
            f"force_torque fell back to a fixed torque centre: '{line}'. The "
            f"zone_id is not linked to an ALE body."
        )
        assert expected in line, (
            f"Expected centre type '{expected}', log says: '{line}'"
        )


def build_case():
    """Load the reference case file and apply the settings common to all tests.

    Callers mutate the returned dict directly for anything test-specific.
    """
    with open(join(HERE, CASE_FILE), "r") as f:
        case = json5.load(f)

    case["case"]["mesh_file"] = join("meshes", MESH_NAME)
    case["case"]["output_directory"] = "./"
    case["case"]["time"]["timestep"] = DT
    case["case"]["time"]["end_time"] = END_TIME
    case["case"]["checkpoint_value"] = CHECKPOINT_TIME
    case["case"].pop("restart_file", None)

    is_dp = conftest.RP == "dp"

    # Solver tolerances.
    solver_tol = 1.0e-13 if is_dp else 1.0e-7
    case["case"]["fluid"]["velocity_solver"]["absolute_tolerance"] = solver_tol
    case["case"]["fluid"]["pressure_solver"]["absolute_tolerance"] = solver_tol

    for comp in case["case"].get("simulation_components", []):
        if comp.get("type") == "force_torque":
            comp["long_print"] = True

    return case


def write_case(case, path):
    with open(path, "w") as f:
        json.dump(case, f, indent=4)
    return path


def recorded_values():
    """Recorded reference for every compared metric, for this precision."""
    return {**REFERENCE_FORCES[conftest.RP], **REFERENCE_KINEMATICS[conftest.RP]}


def analytic_values(t):
    """Closed-form rigid-body kinematics at time `t`."""
    px, py = analytic_pivot_position(t)
    vx, vy = analytic_pivot_velocity(t)
    return {"px": px, "py": py, "vx": vx, "vy": vy,
            "rot_z": analytic_rotation_z_deg(t)}


def log_metrics_table(emit, reference, restarted, step_p1, step_p2):
    """Report every metric next to the values it is checked against.

    The expected column is the recorded reference. 
    Kinematic rows also carry the deviation of the reference run
    from the closed form, which is the discretisation error rather than a
    regression signal.
    """
    t = reference.get(step_p1, {}).get("t")
    expected = recorded_values()
    analytic = analytic_values(t) if t is not None else {}
    width = 118

    def rel(got, want):
        if got is None or want in (None, 0.0):
            return ""
        return f" (rel {abs(got - want) / abs(want):.1e})"

    def cell(name, parsed, step, key, want):
        got = parsed.get(step, {}).get(key)
        if got is None:
            return f"{name} {'MISSING':<19}"
        return f"{name} {got:+.12e}{rel(got, want)}"

    emit("=" * width)
    emit("ALE metrics (%s) at t = %s: reference step %d vs restarted step %d",
         conftest.RP.upper(), f"{t:.6e}" if t is not None else "unknown",
         step_p1, step_p2)
    emit("expected = recorded reference; 'analytic' = part 1 against the "
         "closed form")
    emit("=" * width)
    for key in FORCE_KEYS + KINEMATIC_KEYS:
        want = expected.get(key)
        trailer = ""
        if key in analytic:
            got = reference.get(step_p1, {}).get(key)
            trailer = f" | analytic{rel(got, analytic[key])}"
        emit("[%-8s] expected %-19s | %s | %s%s", key.upper(),
             f"{want:+.12e}" if want is not None else "n/a",
             cell("part1", reference, step_p1, key, want),
             cell("part2", restarted, step_p2, key, want),
             trailer)
    emit("=" * width)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def invocation_dir(request):
    """Directory pytest was started in. Unaffected by the tests' chdir."""
    return Path(request.config.invocation_params.dir)


@pytest.fixture(scope="session")
def ale_neko(tmp_path_factory, invocation_dir):
    """Compile the ALE user file once per session, return the binary path."""
    makeneko = resolve_executable(get_makeneko(), invocation_dir)
    build_dir = tmp_path_factory.mktemp("ale_build")
    shutil.copy(join(HERE, USER_FILE), build_dir)

    result = subprocess.run(
        [makeneko, USER_FILE],
        cwd=build_dir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(
            f"makeneko failed on {USER_FILE} (exit {result.returncode}):\n"
            f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
        )

    binary = join(build_dir, "neko")
    assert isfile(binary), f"makeneko succeeded but produced no binary in {build_dir}"
    return abspath(binary)


@pytest.fixture
def logs_dir(request):
    """Absolute path to the `logs` directory, created if necessary.

    Absolute because the tests chdir into a temporary working directory, and
    because CI collects this directory as an artifact.
    """
    directory = request.config.rootpath / "logs"
    directory.mkdir(exist_ok=True)
    return directory


@pytest.fixture
def ale_workdir(tmp_path, monkeypatch):
    """A clean working directory containing the mesh, and chdir'ed into.

    `monkeypatch.chdir` is restored by pytest even if the test raises.
    """
    os.makedirs(tmp_path / "meshes", exist_ok=True)
    shutil.copy(find_mesh(MESH_NAME), tmp_path / "meshes" / MESH_NAME)
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.fixture
def diag(request):
    """Return `emit(msg, *args)` for the diagnostic tables.

    To see them on a passing run, either:

        pytest --log-cli-level=DEBUG ...    built in, works as-is
        pytest --show-diagnostics ...       needs the option registered in
                                            tests/integration/conftest.py

    The second is optional; if the option is not registered this falls back to
    the default quietly.
    """
    try:
        verbose = bool(request.config.getoption("--show-diagnostics"))
    except ValueError:
        verbose = False  # option not registered in conftest.py
    level = logging.INFO if verbose else logging.DEBUG
    return lambda msg, *args: logger.log(level, msg, *args)


@pytest.fixture
def neko_runner(launcher_script, ale_neko, logs_dir, invocation_dir, request):
    """Return `run(case_file, log_name)`, which fails the test on a non-zero exit."""
    launcher = resolve_executable(launcher_script, invocation_dir)
    assert os.access(launcher, os.X_OK), (
        f"The launcher script {launcher} is not executable. Run "
        f"`chmod +x {launcher}`."
    )

    def run(case_file, log_name, nprocs=2):
        log_path = str(logs_dir / f"{request.node.name}_{log_name}.log")
        result = run_neko(launcher, configure_nprocs(nprocs), case_file,
                          ale_neko, log_path)
        if result.returncode != 0:
            pytest.fail(
                f"neko exited with {result.returncode} on {case_file}.\n"
                f"Last lines of {log_path}:\n{tail(log_path)}"
            )
        return log_path

    return run


# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------
def test_ale(ale_workdir, neko_runner, diag):
    """Rigid-body ALE motion, the user kinematics hook, and checkpoint restart.
    
    """
    run = neko_runner
    tol = TOLERANCES[conftest.RP]

    # -- part 1: reference run, 0 -> END_TIME -------------------------------
    case = build_case()
    write_case(case, "part1.case")
    log_p1 = run("part1.case", "part1")

    # Guard against silently testing a fixed torque centre.
    assert_ale_link_active(log_p1)

    # -- pick the checkpoint written at CHECKPOINT_TIME ---------------------
    checkpoints = sorted(f for f in os.listdir(".") if f.endswith(".chkp"))
    assert len(checkpoints) >= 2, (
        f"Expected at least two .chkp files after part 1, found {checkpoints}. "
        f"Check output_checkpoints/checkpoint_control in the case file."
    )
    restart_file = checkpoints[-2]
    diag("Restarting from %s (of %s)", restart_file, checkpoints)

    # -- part 2: restarted run, CHECKPOINT_TIME -> END_TIME ------------------
    case["case"]["restart_file"] = restart_file
    use_phi_from_first_run(case["case"])
    write_case(case, "part2.case")
    log_p2 = run("part2.case", "part2")

    # -- parse -------------------------------------------------------------
    reference = parse_ale_log(log_p1)
    restarted = parse_ale_log(log_p2)
    assert reference, f"No metrics parsed from {log_p1}:\n{tail(log_p1)}"
    assert restarted, f"No metrics parsed from {log_p2}:\n{tail(log_p2)}"

    # No residual comparison here.
    overlap = overlap_steps(reference, restarted)

    step_p1, step_p2 = COMPARE_STEP_P1, COMPARE_STEP_P2
    assert step_p2 in overlap, (
        f"t={COMPARE_TIME} (reference step {step_p1}, restarted step {step_p2}) "
        f"is not in the overlap window {overlap}."
    )
    log_metrics_table(diag, reference, restarted, step_p1, step_p2)

    require_metrics(reference, step_p1, FORCE_KEYS + KINEMATIC_KEYS, log_p1,
                    "reference run")
    require_metrics(restarted, step_p2, FORCE_KEYS + KINEMATIC_KEYS, log_p2,
                    "restarted run")

    # -- 1. the restart resumed at the right time --------------------------
    for step in overlap:
        assert_allclose(
            restarted[step]["t"], reference[step + RESTART_OFFSET]["t"],
            rtol=0.0, atol=0.1 * DT,
            err_msg=(
                f"Restarted step {step} is at t={restarted[step]['t']}, but "
                f"reference step {step + RESTART_OFFSET} is at "
                f"t={reference[step + RESTART_OFFSET]['t']}. The restart most "
                f"likely resumed from the wrong checkpoint ({restart_file})."
            ),
        )

    # -- 2. restart consistency over the whole overlap window ---------------
    for step in overlap:
        got, want = restarted[step], reference[step + RESTART_OFFSET]
        for key in sorted(set(got) & set(want) - {"t"}):
            assert_allclose(
                got[key], want[key],
                rtol=tol["restart_rtol"], atol=tol["restart_atol"],
                err_msg=(
                    f"Restart mismatch in '{key}': restarted step {step} vs "
                    f"reference step {step + RESTART_OFFSET}"
                ),
            )

    # -- 3. rigid-body kinematics against the closed form -------------------
    for label, parsed, step in (("reference", reference, step_p1),
                                ("restarted", restarted, step_p2)):
        t = parsed[step]["t"]
        px, py = analytic_pivot_position(t)
        vx, vy = analytic_pivot_velocity(t)

        for key, expected in (("px", px), ("py", py), ("vx", vx), ("vy", vy)):
            assert_allclose(
                parsed[step][key], expected,
                rtol=tol["kin_rtol"], atol=tol["kin_atol"],
                err_msg=f"{label} run: {key} at t={t} disagrees with the "
                        f"analytic rigid-body motion",
            )

        assert_allclose(
            parsed[step]["rot_z"], analytic_rotation_z_deg(t),
            rtol=0.0, atol=tol["angle_atol"],
            err_msg=f"{label} run: rotation angle at t={t} disagrees with the "
                    f"analytic rigid-body motion",
        )

        # No motion was prescribed in z.
        if "pz" in parsed[step]:
            assert_allclose(parsed[step]["pz"], PIVOT_0[2],
                            rtol=0.0, atol=tol["kin_atol"],
                            err_msg=f"{label} run: pivot drifted in z")

    # -- 4. fluid forces against the recorded reference ---------------------
    expected_forces = REFERENCE_FORCES[conftest.RP]
    for label, parsed, step in (("reference", reference, step_p1),
                                ("restarted", restarted, step_p2)):
        for key, expected in expected_forces.items():
            assert_allclose(
                parsed[step][key], expected, rtol=tol["force_rtol"],
                err_msg=f"{label} run: {key} deviates from the recorded value",
            )


    # -- 5. kinematics against the recorded reference -----------------------
    expected_kinematics = REFERENCE_KINEMATICS[conftest.RP]
    for label, parsed, step in (("reference", reference, step_p1),
                                ("restarted", restarted, step_p2)):
        for key, expected in expected_kinematics.items():
            assert_allclose(
                parsed[step][key], expected, rtol=tol["kin_recorded_rtol"],
                err_msg=f"{label} run: {key} deviates from the recorded value",
            )


def test_restart_reproduces_residuals(ale_workdir, neko_runner, diag):
    """Restart completeness, probed through the pressure solver residuals.

    This configuration has its own solution, so REFERENCE_FORCES does not
    apply; the test asserts self-consistency only.
    """
    run = neko_runner
    rtol = EXACT_RESTART_RTOL[conftest.RP]
    case = build_case()
    cfg = case["case"]

    cfg["time"]["end_time"] = FAST_END_TIME
    cfg["checkpoint_value"] = FAST_CHECKPOINT_TIME

    # The only deviation from production that matters to the comparison.
    cfg["fluid"]["velocity_solver"]["projection_space_size"] = 0
    cfg["fluid"]["pressure_solver"]["projection_space_size"] = 0

    tol = FAST_SOLVER_TOL[conftest.RP]
    cfg["fluid"]["velocity_solver"]["absolute_tolerance"] = tol
    cfg["fluid"]["pressure_solver"]["absolute_tolerance"] = tol
    cfg["fluid"]["ale"]["solver"]["absolute_tolerance"] = tol

    write_case(case, "exact1.case")
    log_p1 = run("exact1.case", "exact_part1")

    checkpoints = sorted(f for f in os.listdir(".") if f.endswith(".chkp"))
    assert len(checkpoints) >= 2, (
        f"Expected at least two .chkp files, found {checkpoints}."
    )
    cfg["restart_file"] = checkpoints[-2]
    use_phi_from_first_run(cfg)
    write_case(case, "exact2.case")
    log_p2 = run("exact2.case", "exact_part2")

    reference = parse_ale_log(log_p1)
    restarted = parse_ale_log(log_p2)
    overlap = overlap_steps(reference, restarted, offset=FAST_RESTART_OFFSET)

    res_ref = get_field_residuals(parse_residuals(log_p1), PRESSURE_FIELD, log_p1)
    res_rst = get_field_residuals(parse_residuals(log_p2), PRESSURE_FIELD, log_p2)

    compared = 0
    for step in overlap:
        ref_step = step + FAST_RESTART_OFFSET
        if step not in res_rst or ref_step not in res_ref:
            continue
        got, want = res_rst[step], res_ref[ref_step]
        compared += 1

        diag(
            "step %d/%d | res_start %.9e vs %.9e | res_final %.9e vs %.9e | "
            "iters %d vs %d", step, ref_step, got.res_start, want.res_start,
            got.res_final, want.res_final, got.iters, want.iters,
        )

        assert_allclose(
            got.res_start, want.res_start, rtol=rtol,
            err_msg=(
                f"Pressure start residual differs at restarted step {step} "
                f"(reference step {ref_step}). The state handed to the solver "
                f"after the restart is not the state the reference run had at "
                f"the same instant -- something the checkpoint should carry "
                f"was not restored."
            ),
        )
        assert_allclose(
            got.res_final, want.res_final, rtol=rtol,
            err_msg=f"Pressure final residual differs at restarted step {step}",
        )
        # Identical systems should take an identical number of iterations;
        # allow a margin, since the count is a discrete function of the
        # residual and can flip either side of the convergence threshold.
        assert abs(got.iters - want.iters) <= MAX_ITER_DIFF, (
            f"Pressure iteration count differs by more than {MAX_ITER_DIFF} "
            f"at restarted step {step}: {got.iters} vs {want.iters} at "
            f"reference step {ref_step}"
        )

    assert compared >= MIN_OVERLAP, (
        f"Only {compared} steps had residuals in both runs; expected at least "
        f"{MIN_OVERLAP}. Check that '{PRESSURE_FIELD}' is the name the pressure "
        f"solver reports itself under."
    )

    # The fields themselves should agree just as tightly.
    for step in overlap:
        got, want = restarted[step], reference[step + FAST_RESTART_OFFSET]
        for key in sorted(set(got) & set(want) - {"t"}):
            assert_allclose(
                got[key], want[key], rtol=rtol,
                atol=TOLERANCES[conftest.RP]["kin_atol"],
                err_msg=(
                    f"Restart mismatch in '{key}': restarted step {step} vs "
                    f"reference step {step + FAST_RESTART_OFFSET}"
                ),
            )