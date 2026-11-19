import reframe as rfm
import reframe.utility.sanity as sn
import os
import csv
import statistics
import string
import json

def get_gpu_device(partition):
    for device in partition.devices:
        if device.type == 'gpu':
            return device

def set_reference(test, name, ref):
    '''Attach a reference value for `name` on the test's current partition.

    Use this rather than `test.reference.setdefault(partition, {})[name] =
    ref`. `reference` is a ScopedDict: it flattens a nested {scope: {var:
    ref}} mapping into 'scope:var' keys when the whole mapping is *assigned*,
    but mutating the plain dict that setdefault hands back merely stores a
    nested dict under the scope key. check_performance looks up 'scope:var',
    never finds it, falls back to (0, None, None) -- no bounds at all -- and
    the variable then passes unconditionally.

    A gate written the setdefault way therefore appears in the performance
    report and enforces nothing. Confirmed 2026-09-06 by regressing
    field_add2 on purpose: the measured wrapper overhead went to 1.64x
    against a 1.30 cap and ReFrame still reported 'pass'.
    '''
    test.reference[f'{test.current_partition.fullname}:{name}'] = ref

class NekoError(Exception):
    pass

# This class is basically a copy of Autotools from reframe, but with support
# for a configuredir variable.
# Copyright 2016-2022 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# Under BSD-3-Clause license
# Modified by Neko authors 2022
class OutOfSourceAutotools(rfm.core.buildsystems.ConfigureBasedBuildSystem):
    '''A build system for compiling Autotools-based projects.

    This build system will emit the following commands:

    1. Create a build directory if :attr:`builddir` is not :class:`None` and
       change to it.
    2. Invoke ``configure`` to configure the project by setting the
       corresponding flags for compilers and compiler flags.
    3. Issue ``make`` to compile the code.

    ReFrame 4.8 added :attr:`configuredir` to ConfigureBasedBuildSystem
    itself, so redefining it here is now a hard load error -- ReFrame skips
    the whole file with 'variable configuredir is already defined', reports
    'Found 0 check(s)' and still exits 0, which is why this went unnoticed:
    the CI job stayed green while running nothing at all.

    That also makes this class largely redundant: the stock Autotools build
    system supports out-of-source builds since 4.8. Swapping to it would
    change behaviour for every test here, including on Cray systems that
    cannot be exercised from this workspace, so it is left as a follow-up
    rather than folded into this fix.
    '''

    def emit_build_commands(self, environ):
        prepare_cmd = []
        if self.srcdir:
            prepare_cmd += ['cd %s' % self.srcdir]

        if self.builddir:
            prepare_cmd += ['mkdir -p %s' % self.builddir,
                            'cd %s' % self.builddir]

        if self.builddir:
            configure_cmd = [os.path.join(
                os.path.relpath(self.configuredir, self.builddir), 'configure')]
        else:
            configure_cmd = [os.path.join(self.configuredir, 'configure')]

        cc = self._cc(environ)
        cxx = self._cxx(environ)
        ftn = self._ftn(environ)
        cppflags = self._cppflags(environ)
        cflags   = self._cflags(environ)
        cxxflags = self._cxxflags(environ)
        fflags   = self._fflags(environ)
        ldflags  = self._ldflags(environ)
        if cc:
            configure_cmd += ['CC="%s"' % cc]

        if cxx:
            configure_cmd += ['CXX="%s"' % cxx]

        if ftn:
            configure_cmd += ['FC="%s"' % ftn]

        if cppflags:
            configure_cmd += ['CPPFLAGS="%s"' % ' '.join(cppflags)]

        if cflags:
            configure_cmd += ['CFLAGS="%s"' % ' '.join(cflags)]

        if cxxflags:
            configure_cmd += ['CXXFLAGS="%s"' % ' '.join(cxxflags)]

        if fflags:
            configure_cmd += ['FCFLAGS="%s"' % ' '.join(fflags)]

        if ldflags:
            configure_cmd += ['LDFLAGS="%s"' % ' '.join(ldflags)]

        if self.config_opts:
            configure_cmd += self.config_opts

        make_cmd = ['make -j']
        if self.max_concurrency is not None:
            make_cmd += [str(self.max_concurrency)]

        if self.make_opts:
            make_cmd += self.make_opts

        return prepare_cmd + [' '.join(configure_cmd), ' '.join(make_cmd)]

class BuildNeko(rfm.CompileOnlyRegressionTest):
    #build_system = 'Autotools'
    build_system = OutOfSourceAutotools()
    builddir = 'build'
    backend = variable(str)
    real = parameter(os.getenv('NEKO_REAL', 'dp,sp').split(','))

    @run_after('setup')
    def set_backend(self):
        gpu_device = get_gpu_device(self.current_partition)
        if gpu_device is None:
            self.backend = 'cpu'
        else:
            self.backend = 'device'
            self.gpu_device = gpu_device

    @run_before('compile')
    def prepare_build(self):
        self.build_system.configuredir = os.path.join(self.prefix, '../../')

        self.build_system.max_concurrency = 32
        self.build_system.make_opts = ['install']

        self.install_dir = os.path.join(self.stagedir, 'install')
        self.build_system.config_opts.append(f'--prefix={self.install_dir}')

        self.build_system.config_opts.append(f'--enable-real={self.real}')

        if self.backend == 'device':
            config = ''
            if self.gpu_device.arch == 'amd':
                config = '--with-hip="$HIP_PATH"'
            elif self.gpu_device.arch == 'nvidia':
                config = '--with-cuda="$CUDA_HOME"'
                # CUDA_ARCH is substituted verbatim into the nvcc command
                # line (src/Makefile.am), so it has to carry its own flag:
                # '-arch=sm_86', not 'sm_86'. A bare value makes nvcc fail
                # with 'A single input file is required for a non-link
                # phase'. The compute capability is a property of the
                # machine, so it lives in the partition's extras next to
                # select_device rather than being hardcoded here.
                cuda_arch = self.current_partition.extras.get('cuda_arch')
                if cuda_arch:
                    self.build_system.config_opts.append(
                        f'CUDA_ARCH="{cuda_arch}"')
            else:
                raise NekoError(f'Unknown gpu arch {self.gpu_device.arch}')

            self.build_system.config_opts.append(config)

    @sanity_function
    def validate_build(self):
        config = os.path.join(self.stagedir, 'src', 'config', 'neko_config.f90')
        if self.backend == 'cpu':
            return sn.assert_not_found(r'NEKO_BCKND_\w+ = 1', config)
        else:
            return sn.assert_found(r'NEKO_BCKND_\w+ = 1', config)

# Use this for children of NekoTestBase that don't need makeneko
class DummyBuildSystem(rfm.core.buildsystems.BuildSystem):
    def emit_build_commands(self, environ):
        return []

class MakeNeko(rfm.core.buildsystems.BuildSystem):
    srcfile = variable(str, type(None), value=None)

    def __init__(self, neko_build):
        self.makeneko = os.path.join(neko_build.install_dir, 'bin', 'makeneko')

    def emit_build_commands(self, environ):
        if not self.srcfile:
            raise NekoError('Source file required')

        return [f'{self.makeneko} "{self.srcfile}"']

class NekoTestBase(rfm.RegressionTest):
    valid_systems = ['*']
    valid_prog_environs = ['PrgEnv-cray', 'PrgEnv-gnu', 'PrgEnv-intel','default']
    neko_build = fixture(BuildNeko, scope='environment')

    scheme = parameter(os.getenv('NEKO_SCHEME', 'pnpn').split(','))
    case = variable(str)

    mesh_file = variable(str, value='')
    dt = variable(float, value=0.0)
    T_end = variable(float, value=0.0)

    abstol_vel = {'sp': 1e-5, 'dp': 1e-9}
    abstol_prs = {'sp': 1e-5, 'dp': 1e-9}

    # Set dofs to enable workrate perf var
    dofs = variable(int, value=0)
    first_workrate_timestep = variable(int, value=0)

    @run_before('compile')
    def copy_mesh_file(self):
        if self.mesh_file == '':
            return

        src = os.path.join(self.prefix, '..', self.mesh_file)
        dst = os.path.join(self.stagedir, self.mesh_file)
        self.postbuild_cmds += [
                f'mkdir -p {os.path.dirname(self.mesh_file)}',
                f'cp "{src}" "{dst}"'
        ]

    @run_before('run')
    def make_case_file(self):
        case_file = os.path.join(self.stagedir, self.case)
        case_template = case_file + '.template'

        self.executable_opts.append(self.case)

        if os.path.exists(case_file):
            pass
        elif os.path.exists(case_template):
            with open(case_template) as tf:
                case_json = json.load(tf)
            case_json["case"]["fluid"]["velocity_solver"]["absolute_tolerance"] = \
                self.abstol_vel[self.neko_build.real]
            case_json["case"]["fluid"]["pressure_solver"]["absolute_tolerance"] = \
                self.abstol_prs[self.neko_build.real]
            case_json["case"]["fluid"]["scheme"] = self.scheme
            case_json["case"]["mesh_file"] = self.mesh_file
            case_json["case"]["timestep"] = self.dt
            case_json["case"]["end_time"] = self.T_end

            with open(case_file, 'w') as cf:
                json.dump(case_json, cf, indent=2)
        else:
            raise NekoError(f'Cannot find {case_file} or {case_template}')

    @run_before('run')
    def set_num_tasks(self):
        if self.neko_build.backend == 'cpu':
            num_cpus = self.current_partition.processor.num_cpus
            cpus_per_core = self.current_partition.processor.num_cpus_per_core
            self.num_tasks = int(num_cpus / cpus_per_core)
        elif self.neko_build.backend == 'device':
            gpu_device = get_gpu_device(self.current_partition)
            if gpu_device is None:
                raise NekoError("Device of type gpu not defined for partition!")
            self.num_tasks = gpu_device.num_devices
        else:
            raise NekoError(f'Unknown backend {self.neko_build.backend}!')

    @run_before('run')
    def select_device(self):
        try:
            select_device = self.current_partition.extras['select_device']
            self.executable_opts.insert(0, self.executable)
            self.executable = select_device
        except KeyError:
            pass

    @sanity_function
    def normal_end(self):
        return sn.assert_found('Normal end.', self.stdout)

    @run_before('performance')
    def set_time_perf(self):
        timesteps = sn.extractall(r'Elapsed time \(s\):\s+(\S+)', self.stdout, 1, float)

        pf = sn.make_performance_function(lambda: timesteps[-1], 's')
        self.perf_variables['total_runtime'] = pf

        if self.dofs != 0:
            pes = self.num_tasks

            def workrate():
                end = sn.count(timesteps) - 1
                time = timesteps[end] - timesteps[self.first_workrate_timestep]
                dofs = 8**3 * 32**3
                iters = end - self.first_workrate_timestep
                return 1e-3 * dofs * iters / time / pes

            pf = sn.make_performance_function(workrate, 'Mdofs/s/pe')
            self.perf_variables['workrate'] = pf

class GetTgvDns(rfm.RunOnlyRegressionTest):
    descr = 'Download TGV DNS data'
    executable = './get-tgv-dns.sh'
    local = True

    @run_after('run')
    def load_enstrophy(self):
        self.enstrophy = {}
        path = os.path.join(self.stagedir, 'spectral_Re1600_512.gdiag')
        with open(path, newline='') as f:
            reader = csv.reader(f, delimiter=' ')
            for row in reader:
                if row[0][0] == '#':
                    continue
                # time: value
                self.enstrophy[float(row[0])] = float(row[3])

    @sanity_function
    def check_data_count(self):
        return sn.assert_eq(sn.count(sn.defer(self.enstrophy)), 2000)

class TgvBase(NekoTestBase):
    descr = 'Run TGV and compare with DNS data'
    executable = './neko'
    case = 'tgv.case'
    tgv_dns = fixture(GetTgvDns, scope='session')

    @run_after('setup')
    def set_build(self):
        self.build_system = MakeNeko(self.neko_build)
        self.sourcepath = 'tgv.f90'

    @sn.deferrable
    def max_error(self, time_ens):
        errs = []
        for time, ens in time_ens:
            # Round time to 3 decimals to find corresponding DNS sample
            time = round(time, 3)
            if time == 20.0:
                # DNS data does not include the last timestep
                continue
            try:
                dns = self.tgv_dns.enstrophy[time]
            except KeyError:
                raise NekoError(f'DNS enstrophy not sampled at {time}')
            errs.append(100 * abs(1 - ens/dns))
        return max(errs)

    @performance_function('%')
    def enstrophy_error(self):
        time_ens = sn.extractall(r'Time: (\S+).*Enstrophy: (\S+)', self.stdout, (1, 2), (float, float))
        return self.max_error(time_ens)

@rfm.simple_test
class Tgv8(TgvBase):
    mesh_file = '../examples/tgv/512.nmsh'
    dt = 1e-2
    T_end = 20.0

    @run_before('performance')
    def set_reference(self):
        if self.neko_build.real == 'dp':
            self.reference = {
                'dt:gpu': {
                    'total_runtime': (45, -0.50, 0.10, 's'),
                },
                'dt:cpu': {
                    'total_runtime': (16, -0.50, 0.10, 's'),
                },
            }

            # For all systems.
            #
            # NOTE: this reference was previously attached with
            # `self.reference.setdefault(...)[...] = ...`, which silently
            # attaches nothing -- see set_reference() above. The accuracy
            # gate has therefore been inert, so this is the first time it is
            # actually enforced. If a run starts failing on
            # enstrophy_error, the value below is the thing to re-measure,
            # not the mechanism.
            set_reference(self, 'enstrophy_error',
                          (33.48, -0.01, 0.01, '%'))

@rfm.simple_test
class Tgv32(TgvBase):
    mesh_file = '../examples/tgv/32768.nmsh'
    dt = 1e-3
    T_end = 20.0
    dofs = 8**3 * 32**3
    # Where flow has become turbulent
    first_workrate_timestep = 12000

    @run_before('performance')
    def set_reference(self):
        if self.neko_build.real == 'dp':
            self.reference = {
                'dt:gpu': {
                    'total_runtime': (4800, -0.50, 0.05, 's'),
                }
            }

            # For all systems. Inert until now for the same reason as Tgv8's
            # -- see the note there and set_reference() above.
            set_reference(self, 'enstrophy_error', (6.73, -0.01, 0.01, '%'))

@rfm.simple_test
class MiniHemi(NekoTestBase):
    descr = 'Two iterations of hemi as a smoke test'
    build_system = DummyBuildSystem()
    case = 'minihemi.case'
    mesh_file = '../examples/hemi/hemi.nmsh'

    @run_before('compile')
    def setup_case(self):
        self.executable = os.path.join(self.neko_build.install_dir, 'bin/neko')

@rfm.simple_test
class MiniTgv8(NekoTestBase):
    descr = 'Two iterations of TGV as a smoke test'
    mesh_file = '../examples/tgv/512.nmsh'
    dt = 1e-2
    T_end = 0.02
    executable = './neko'
    case = 'tgv.case'

    @run_after('setup')
    def set_build(self):
        self.build_system = MakeNeko(self.neko_build)
        self.sourcepath = 'tgv.f90'

@rfm.simple_test
class MiniRB(NekoTestBase):
    descr = 'Two iterations of 3D RB as a smoke test'
    mesh_file = '../examples/rayleigh_benard/box.nmsh'
    dt = 1e-2
    T_end = 0.02
    executable = './neko'
    case = 'rayleigh.case'

    @run_after('setup')
    def set_build(self):
        self.build_system = MakeNeko(self.neko_build)
        self.sourcepath = 'rayleigh.f90'

    # Restrict small case to 2 tasks
    @run_before('run')
    def set_num_tasks(self):
        if self.neko_build.backend == 'cpu':
            self.num_tasks = 2

# ---------------------------------------------------------------------------
# math_ops: throughput of math/device_math, and the cost of the wrappers
# ---------------------------------------------------------------------------

@sn.deferrable
def _fmean(values):
    return statistics.fmean(sn.evaluate(values))

@sn.deferrable
def _fstdev(values):
    return statistics.stdev(sn.evaluate(values))

class MakeBench(rfm.core.buildsystems.BuildSystem):
    '''Build a tests/bench driver against the fixture's Neko install.

    MakeNeko cannot be used for these: makeneko refuses any input without a
    `module user` or a type-injecting module, and it generates its own
    program, whereas the bench drivers are standalone programs of their own.
    The committed bench Makefile resolves the compiler, flags and libraries
    through pkg-config, so pointing PKG_CONFIG_PATH at the install is the
    whole of the integration.

    `make clean` runs first because ReFrame copies the source directory
    verbatim into the stage dir. Build products there are gitignored, so a
    developer who has built the benchmark by hand would otherwise ship stale
    objects into the test and never see it.
    '''

    def __init__(self, neko_build):
        self.pkgconfig = os.path.join(neko_build.install_dir, 'lib',
                                      'pkgconfig')

    def emit_build_commands(self, environ):
        return [
            'make clean',
            f'PKG_CONFIG_PATH={self.pkgconfig}:$PKG_CONFIG_PATH make'
        ]

class MathOpsBase(rfm.RegressionTest):
    '''Shared setup for the tests/bench/math_ops benchmark.

    Does not derive from NekoTestBase: that class builds a case file from a
    template and asserts 'Normal end.', which simulation.f90 prints and this
    driver never does. Only the BuildNeko fixture is shared.
    '''

    valid_systems = ['*']
    valid_prog_environs = ['PrgEnv-cray', 'PrgEnv-gnu', 'PrgEnv-intel',
                           'default']
    neko_build = fixture(BuildNeko, scope='environment')

    sourcesdir = os.path.join('..', 'bench', 'math_ops')
    executable = './mathbench'

    #: Fixed deliberately -- the glsc3 references below are values *of this
    #: mesh*, so changing it invalidates them. Relative to tests/.
    mesh_file = 'bench/nekbone/data/512.nmsh'

    #: Must match lx_sweep in tests/bench/math_ops/driver.f90.
    lx_sweep = [2, 3, 4, 6, 8, 12]
    ops = ['add2', 'col2', 'glsc3']
    wrappers = ['field_math', 'vector_math', 'matrix_math']

    #: Iterations per timed loop; 0 means verify only.
    niter = variable(int, value=0)

    @run_after('setup')
    def set_build(self):
        self.build_system = MakeBench(self.neko_build)

    @run_before('compile')
    def copy_mesh_file(self):
        src = os.path.join(self.prefix, '..', self.mesh_file)
        self.postbuild_cmds += [f'cp "{src}" mesh.nmsh']

    @run_before('run')
    def set_executable_opts(self):
        self.executable_opts = ['mesh.nmsh', str(self.niter)]

    @run_before('run')
    def set_runtime_libs(self):
        # Neko links json-fortran as a shared library, so the run needs its
        # directory on LD_LIBRARY_PATH even though the build resolved it
        # through pkg-config. Derived from the .pc file rather than
        # hardcoded, since the dependency prefixes differ per machine.
        pkgconfig = os.path.join(self.neko_build.install_dir, 'lib',
                                 'pkgconfig')
        self.prerun_cmds += [
            f'export PKG_CONFIG_PATH={pkgconfig}:$PKG_CONFIG_PATH',
            'export LD_LIBRARY_PATH='
            '$(pkg-config --libs-only-L neko | sed -e "s/-L//g" '
            '-e "s/ \\+/:/g"):$LD_LIBRARY_PATH'
        ]

    @sanity_function
    def validate_run(self):
        # The driver aborts on any cross-path disagreement, so reaching the
        # end of the sweep is itself the correctness result. The backend
        # assertion catches a device build that silently ran on the host --
        # that would satisfy every numerical check and prove nothing.
        expected_bcknd = 1 if self.neko_build.backend == 'device' else 0
        return sn.all([
            sn.assert_eq(sn.count(sn.findall(r'# verify OK', self.stdout)),
                         len(self.lx_sweep)),
            sn.assert_found(rf'# bcknd_dev : {expected_bcknd}', self.stdout)
        ])

@rfm.simple_test
class MathOpsVerify(MathOpsBase):
    descr = 'math_ops cross-path and rank-count verification'
    niter = 0
    nranks = parameter([1, 2, 4])

    #: Reduced glsc3 by lx, for mesh_file at dp, measured 2026-09-06.
    #:
    #: The fill is a function of the dofmap's global coordinates, so the
    #: global multiset of values -- and hence this reduction -- does not
    #: depend on how the mesh is partitioned. Pinning the values therefore
    #: gates rank-count invariance, which no single run can check against
    #: itself. Observed spread across -np 1/2/4/8 is <= 2.4e-14 relative,
    #: from floating-point reduction not being associative; the 1e-12
    #: tolerance below leaves roughly 40x headroom for a different compiler
    #: or vectorisation.
    #:
    #: Regenerate after any change to the mesh, the fill, or the sweep:
    #:   cd tests/bench/math_ops && make
    #:   mpirun -np 1 ./mathbench ../nekbone/data/512.nmsh 0 | grep GLSC3
    glsc3_ref = {
        2: 0.3464241394253017e+06,
        3: 0.1168359770472830e+07,
        4: 0.2769006114889988e+07,
        6: 0.9344507979646834e+07,
        8: 0.2214916174351548e+08,
        12: 0.7475118853765911e+08,
    }
    glsc3_tol = 1e-12

    @run_before('run')
    def set_num_tasks(self):
        self.num_tasks = self.nranks
        # OpenMPI refuses to start more ranks than cores; srun has no such
        # flag, so this must not be applied blindly.
        launcher = getattr(self.job.launcher, 'registered_name', None)
        processor = self.current_partition.processor
        ncpus = getattr(processor, 'num_cpus', None)
        if launcher == 'mpirun' and ncpus and self.num_tasks > ncpus:
            self.job.launcher.options += ['--oversubscribe']

    @run_before('performance')
    def set_glsc3_perf(self):
        for lx in self.lx_sweep:
            patt = rf'GLSC3 lx={lx} value=\s*(\S+)'
            name = f'glsc3_lx{lx}'
            self.perf_variables[name] = sn.make_performance_function(
                sn.extractsingle(patt, self.stdout, 1, float), ''
            )
            # sp reproduces a different number entirely, and to far fewer
            # digits -- logged for history, not gated.
            if self.neko_build.real == 'dp':
                set_reference(self, name, (self.glsc3_ref[lx],
                                           -self.glsc3_tol,
                                           self.glsc3_tol, ''))

@rfm.simple_test
class MathOpsPerf(MathOpsBase):
    '''Throughput of math/device_math, and what the wrappers cost on top.

    Single-rank on purpose: the question is dispatch cost, and glsc3's
    Allreduce would otherwise dominate it. Rank-count behaviour is
    MathOpsVerify's job.
    '''

    descr = 'math_ops wrapper dispatch overhead and math workrate'
    niter = 200
    num_tasks = 1

    #: Cap on the mean of the per-lx ratio t_wrapper/t_math, as an excess
    #: over parity: 0.30 means the wrappers may not average more than 1.30x
    #: the direct call. Calibrated 2026-09-06 on a -g -O2 build: 8 idle runs
    #: peaked at 1.033, and 6 runs pinned to 2 cores against 2 competing
    #: busy loops peaked at 1.117. The gate is deliberately loose relative to
    #: those, because the regressions it can actually resolve are gross ones
    #: -- an added copy or allocation in a wrapper roughly doubles the ratio,
    #: while a lost inline is a few nanoseconds against a 8 us call and is
    #: invisible at any threshold. Preferring a loose gate that never flakes
    #: over a tight one that cries wolf is the deliberate trade.
    overhead_mean_max = 0.30

    #: Absolute cap on the stddev of the per-lx ratio. Same runs: 0.062 idle,
    #: 0.218 contended. This catches a regression confined to a single size,
    #: which widens the spread while barely moving the mean.
    overhead_spread_cap = 0.35

    def ratios(self, op, wrapper):
        patt = (rf'RATIO op={op} path={wrapper} lx=\d+ value=\s*(\S+)'
                r' value_mean=')
        return sn.extractall(patt, self.stdout, 1, float)

    @sanity_function
    def validate_run(self):
        # Assert the record counts explicitly, because nothing downstream
        # will. ReFrame catches any exception from evaluating a performance
        # variable, logs 'skipping evaluation of performance variable' at
        # warning level and *continues* (pipeline.py, check_performance), so
        # a regex that stops matching does not fail the test -- the metric
        # just silently disappears from the report and the gate stops
        # existing. This was not hypothetical: it swallowed a broken
        # workrate pattern here during development, and the run went green.
        nbench = len(self.ops) * (len(self.wrappers) + 1) * len(self.lx_sweep)
        nratios = len(self.ops) * len(self.wrappers) * len(self.lx_sweep)
        return sn.all([
            super().validate_run(),
            sn.assert_eq(sn.count(sn.findall(r'^BENCH ', self.stdout)),
                         nbench),
            sn.assert_eq(sn.count(sn.findall(r'^RATIO ', self.stdout)),
                         nratios)
        ])

    @run_before('performance')
    def set_overhead_perf(self):
        for op in self.ops:
            for wrapper in self.wrappers:
                r = self.ratios(op, wrapper)
                mean_name = f'overhead_mean_{op}_{wrapper}'
                spread_name = f'overhead_spread_{op}_{wrapper}'
                self.perf_variables[mean_name] = \
                    sn.make_performance_function(_fmean(r), '')
                self.perf_variables[spread_name] = \
                    sn.make_performance_function(_fstdev(r), '')
                # Ratios cancel machine speed, compiler and precision, so
                # unlike the absolute workrate below these are gated
                # everywhere, including on shared CI runners.
                #
                # ReFrame thresholds are relative to the reference, so
                # (1.0, None, 0.30) means 'at most 1.30x, no lower bound' --
                # a wrapper coming out faster is not a regression. The
                # spread has no natural reference value, so it is expressed
                # as an absolute cap: (cap, None, 0.0) means 'at most cap'.
                set_reference(self, mean_name,
                              (1.0, None, self.overhead_mean_max, ''))
                set_reference(self, spread_name,
                              (self.overhead_spread_cap, None, 0.0, ''))

    @run_before('performance')
    def set_workrate_perf(self):
        # The driver always labels the direct path 'math'; which module
        # actually ran is fixed by the build, so name the variable after the
        # backend to keep a perflog row self-describing.
        backend = ('device_math' if self.neko_build.backend == 'device'
                   else 'math')
        for op in self.ops:
            for lx in self.lx_sweep:
                # Every \s* is load-bearing: the driver writes these with an
                # e17.10 edit descriptor, which left-pads a positive value
                # with a blank, so 'min=\S+' does not match 'min= 0.43E-03'.
                patt = (rf'BENCH op={op} path=math lx={lx} n=\d+ '
                        r'min=\s*\S+ mean=\s*\S+ sd=\s*\S+ mdofs=\s*(\S+)')
                self.perf_variables[f'workrate_{backend}_{op}_lx{lx}'] = \
                    sn.make_performance_function(
                        sn.extractsingle(patt, self.stdout, 1, float),
                        'Mdofs/s/pe')

        # No references are set for these. Absolute throughput is a property
        # of the machine, so a reference is only meaningful on a dedicated
        # one -- Dardel, as with Tgv8/Tgv32's total_runtime. Those numbers
        # have to be measured there; they are deliberately not guessed here.
        # See tests/reframe/README.md for how to fill them in.
