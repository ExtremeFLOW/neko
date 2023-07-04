import reframe as rfm
import reframe.utility.sanity as sn
import os
import csv
import string
import json

def get_gpu_device(partition):
    for device in partition.devices:
        if device.type == 'gpu':
            return device

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
    '''

    configuredir = variable(str, value='.')

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
        self.build_system.configuredir = os.path.join(self.prefix, '..')

        self.build_system.max_concurrency = 32
        self.build_system.make_opts = ['install']

        self.install_dir = os.path.join(self.stagedir, 'install')
        self.build_system.config_opts.append(f'--prefix={self.install_dir}')

        self.build_system.config_opts.append(f'--enable-real={self.real}')

        if self.backend == 'device':
            config = ''
            if self.gpu_device.arch == 'amd':
                config = '--with-hip="$HIP_PATH"'
            else:
                raise NekoError(f'Unknown gpu arch {self.gpu_device.arch}')

            self.build_system.config_opts.append(config)

    @sanity_function
    def validate_build(self):
        config = os.path.join(self.stagedir, 'src', 'config', 'neko_config.f90')
        if self.backend == 'cpu':
            return sn.assert_not_found('NEKO_BCKND_\w+ = 1', config)
        else:
            return sn.assert_found('NEKO_BCKND_\w+ = 1', config)

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
    dt = variable(float, value=0)
    T_end = variable(float, value=0)

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
                case_json.dump(cf)
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
    mesh_file = 'examples/tgv/512.nmsh'
    dt = 1e-2
    T_end = 20
    
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

            # For all systems
            ref = self.reference.setdefault(self.current_partition.fullname, {})
            ref['enstrophy_error'] = (33.48, -0.01, 0.01, '%')

@rfm.simple_test
class Tgv32(TgvBase):
    mesh_file = 'examples/tgv/32768.nmsh'
    dt = 1e-3
    T_end = 20    
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

            # For all systems
            ref = self.reference.setdefault(self.current_partition.fullname, {})
            ref['enstrophy_error'] = (6.73, -0.01, 0.01, '%')

@rfm.simple_test
class MiniHemi(NekoTestBase):
    descr = 'Two iterations of hemi as a smoke test'
    build_system = DummyBuildSystem()
    case = 'minihemi.case'
    mesh_file = 'examples/hemi/hemi.nmsh'

    @run_before('compile')
    def setup_case(self):
        self.executable = os.path.join(self.neko_build.install_dir, 'bin/neko')

@rfm.simple_test
class MiniTgv8(NekoTestBase):
    descr = 'Two iterations of TGV as a smoke test'
    mesh_file = 'examples/tgv/512.nmsh'
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
    mesh_file = 'examples/rayleigh-benard/box.nmsh'
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
