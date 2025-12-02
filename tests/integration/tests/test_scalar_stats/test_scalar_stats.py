from os.path import join
from testlib import get_neko, run_neko, configure_nprocs, get_makeneko
import subprocess
import numpy as np
from pysemtools.datatypes.field import FieldRegistry
from pysemtools.datatypes.msh import Mesh
from pysemtools.io.ppymech.neksuite import pynekread
from mpi4py import MPI


def test_scalar_stats(launcher_script, request, log_file, tmp_path):
    """
    Synthetic test for basic scalar statistics.

    This module defines simple bounded, periodic fields u(t), v(t), w(t), and
    p(t) for use in verifying the correctness of a statistics. The fields are
    constructed so that all first- and second-order raw moments (<s>, <us>,
    <vs>, <ws>, <ss>) have known, closed-form analytical values. Note that that
    the fields have no spatial variation, so they are not suitable for testing
    the full statistics output, only the basic one.

    All fields use a common base frequency ω and are defined as

        u(t) = U0 + U1 * cos(ω t)
        v(t) = V0 + V1 * cos(ω t + φ)
        w(t) = W0 + W1 * cos(ω t + ψ)
        s(t) = S0 + S1 * sin(ω t + χ)

    Because each signal is purely harmonic, averaging over an integer number of
    periods yields analytical time-averaged statistics:

    First-order statistics:
        <s> = S0
    Second-order diagonal moments:
        <ss> = S0**2 + 0.5 * S1**2
    Cross moments:
        <us> = U0*S0 + 0.5 * U1*S1 * sin(χ)
        <vs> = V0*S0 + 0.5 * V1*S1 * sin(χ - φ)
        <ws> = W0*S0 + 0.5 * W1*S1 * sin(χ - ψ)


    The test runs 3 scalar_stats simcomps one with xy averaging and one with just
    x. All outputs are tested against analytical values.
    """

    # Gets the nameof the test, i.e. test_demo here. `request` can be use for
    # other things like this.
    test_name = request.node.name

    # Get the path to the neko executable
    neko = "./neko"
    makeneko = get_makeneko()

    result = subprocess.run(
        [makeneko, join("tests", "test_scalar_stats", "test_scalar_stats.f90")],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True)

    assert (
        result.returncode == 0
    ), f"makeneko process failed with exit code {result.returncode}"

    result = subprocess.run(
        ["genmeshbox", "0", "1", "0", "1", "0", "1", "3", "3", "3",
         ".true.", ".true.", ".true."],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True)
    assert (
        result.returncode == 0
    ), f"genmeshbox process failed with exit code {result.returncode}"

    # Number of ranks to launch on
    max_nprocs = 1

    nprocs = configure_nprocs(max_nprocs)

    case_file = join("tests", "test_scalar_stats", "test_scalar_stats.case")

    # Run Neko
    result = run_neko(launcher_script, nprocs, case_file, neko, log_file)

    # Check if the process completed successfully
    assert (
        result.returncode == 0
    ), f"neko process failed with exit code {result.returncode}"

    # Analytical correct values
    U0 = 1.0
    U1 = 1.0
    V0 = 2.0
    V1 = 2.0
    W0 = 3.0
    W1 = 3.0
    S0 = 4.0
    S1 = 4.0
    omega = 4.0
    phi = 3.14192 / 4.0
    psi = 3.14192 / 2.0
    xi = 3.14192 / 3.0

    correct = np.array([
        S0,
        U0 * S0 + 0.5 * U1 * S1 * np.sin(xi),
        V0 * S0 + 0.5 * V1 * S1 * np.sin(xi - phi),
        W0 * S0 + 0.5 * W1 * S1 * np.sin(xi - psi),
        S0**2 + 0.5 * S1**2,
    ])

    #
    # 1D statistics output
    #


    csv = np.genfromtxt(join("tests", "test_scalar_stats", "scalar_stats0.csv"),
                        delimiter=",")[12, 2:]

    quants = ["<s>", "<us>", "<vs>", "<ws>", "<s^2>"]

    for i, q in enumerate(quants):
        error = (csv[i] - correct[i]) / correct[i]
        assert (
            error < 1e-2
        ), f"Error in {q}, {csv[i]} exceeded tolerance: {error}"


    #
    # 2D statistics output
    #
    comm = MPI.COMM_WORLD

    mesh = Mesh(comm, create_connectivity=True)
    fld = FieldRegistry(comm)
    pynekread(join("tests", "test_scalar_stats", "stats2d0.f00001"), comm,
              data_dtype=np.single, msh=mesh, fld=fld)

    # Ordering of keys corresponds to 'correct' array above
    ordered_keys = ["p", "u", "v", "s0", "t"]

    # We check at an essentially random node, the should all have the same
    # values
    for i, key in enumerate(ordered_keys):
        val = fld.registry[key][1,0,0,0]
        assert(
            (val - correct[i]) / correct[i] < 1e-2
        ), f"Value for {key} is {val}, should be {correct[i]}, \
             error exceeded tolerance: {(val - correct[i]) / correct[i]}"

    #
    # 3D statistics output
    #

    mesh = Mesh(comm, create_connectivity=True)
    fld = FieldRegistry(comm)
    pynekread(join("tests", "test_scalar_stats", "stats3d0.f00001"), comm,
              data_dtype=np.single, msh=mesh, fld=fld)

    # Ordering of keys corresponds to 'correct' array above
    ordered_keys = ["p", "u", "v", "w", "t"]

    # We check at an essentially random node, the should all have the same
    # values
    for i, key in enumerate(ordered_keys):
        val = fld.registry[key][1,0,0,0]
        assert(
            (val - correct[i]) / correct[i] < 1e-2
        ), f"Value for {key} is {val}, should be {correct[i]}, \
             error exceeded tolerance: {(val - correct[i]) / correct[i]}"


