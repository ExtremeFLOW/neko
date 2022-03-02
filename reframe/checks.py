import reframe as rfm
import reframe.utility.sanity as sn
import os
import csv
import string

def get_gpu_device(partition):
    for device in partition.devices:
        if device.type == 'gpu':
            return device

class NekoError(Exception):
    pass

class BuildNeko(rfm.CompileOnlyRegressionTest):
    build_system = 'Autotools'
    backend = variable(str)
    sourcesdir = '..'
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
        self.build_system.max_concurrency = 32
        self.build_system.make_opts = ['install']
        self.install_dir = os.path.join(self.stagedir, 'install')
        self.build_system.config_opts.append(f'--prefix={self.install_dir}')
        self.build_system.config_opts.append(f'--enable-real={self.real}')

        self.prebuild_cmds = [
            './regen.sh'
        ]

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
    valid_systems = ['dt:cpu', 'dt:gpu']
    valid_prog_environs = ['cray']
    neko_build = fixture(BuildNeko, scope='environment')
    sourcesdir = '..'

    scheme = parameter(os.getenv('NEKO_SCHEME', 'plan4,pnpn').split(','))
    case = variable(str)

    mesh_file = variable(str, value='')
    dt = variable(str, value='')

    abstol_vel = {'sp': '1d-5', 'dp': '1d-9'}
    abstol_prs = {'sp': '1d-5', 'dp': '1d-9'}

    # Set dofs to enable workrate perf var
    dofs = variable(int, value=0)
    first_workrate_timestep = variable(int, value=0)

    @run_before('run')
    def make_case_file(self):
        case_file = os.path.join(self.stagedir, self.case)
        case_template = case_file + '.template'

        self.executable_opts.append(self.case)

        if os.path.exists(case_file):
            pass
        elif os.path.exists(case_template):
            with open(case_template) as tf:
                ts = tf.read()
            template = string.Template(ts)

            keys = {
                'abstol_vel': self.abstol_vel[self.neko_build.real],
                'abstol_prs': self.abstol_prs[self.neko_build.real],
                'fluid_scheme': self.scheme,
                'mesh_file': self.mesh_file,
                'dt': self.dt,
            }

            ss = template.substitute(keys)
            with open(case_file, 'w') as cf:
                cf.write(ss)
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
    def dt_select_gpu(self):
        if self.current_partition.fullname == 'dt:gpu':
            self.executable_opts.insert(0, self.executable)
            self.executable = 'reframe/rocm-select-gpu'

    @sanity_function
    def normal_end(self):
        return sn.assert_found('normal end.', self.stdout)

    @run_before('performance')
    def set_time_perf(self):
        timesteps = sn.extractall(r'Elapsed time \(s\):\s+(\S+)', self.stdout, 1, float)

        pf = sn.make_performance_function(lambda: timesteps[-1], 's')
        self.perf_variables['total_runtime'] = pf

        if self.dofs != 0:
            def workrate():
                end = sn.count(timesteps) - 1
                time = timesteps[end] - timesteps[self.first_workrate_timestep]
                dofs = 8**3 * 32**3
                iters = end - self.first_workrate_timestep
                return 1e-3 * dofs * iters / time

            pf = sn.make_performance_function(workrate, 'Mdofs/s')
            self.perf_variables['workrate'] = pf

class TgvDns(rfm.RunOnlyRegressionTest):
    descr = 'Download TGV DNS data'
    sourcesdir = '.'
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

class Tgv(NekoTestBase):
    descr = 'Run TGV and compare with DNS data'
    executable = './neko'
    case = 'reframe/tgv.case'
    sourcepath = 'bench/tgv32/tgv.f90'
    tgv_dns = fixture(TgvDns, scope='session')

    @run_after('setup')
    def set_build(self):
        self.build_system = MakeNeko(self.neko_build)

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
class Tgv8(Tgv):
    mesh_file = 'examples/tgv/512.nmsh'
    dt = '1d-2'

    @run_before('performance')
    def set_reference(self):
        if self.neko_build.real == 'dp':
            self.reference = {
                'dt:gpu': {
                    'total_runtime': (120, -0.05, 0.05, 's'),
                    'enstrophy_error': (33.5, -0.01, 0.01, '%'),
                },
                'dt:cpu': {
                    'total_runtime': (16, -0.05, 0.05, 's'),
                    'enstrophy_error': (33.48, -0.01, 0.01, '%'),
                },
            }

@rfm.simple_test
class Tgv32(Tgv):
    mesh_file = 'examples/tgv/32768.nmsh'
    dt = '1d-3'
    dofs = 8**3 * 32**3
    # Where flow has become turbulent
    first_workrate_timestep = 12000

    @run_before('performance')
    def set_reference(self):
        if self.neko_build.real == 'dp':
            self.reference = {
                'dt:gpu': {
                    'total_runtime': (4980, -0.05, 0.05, 's'),
                    'enstrophy_error': (6.73, -0.01, 0.01, '%'),
                }
            }

@rfm.simple_test
class MiniHemi(NekoTestBase):
    descr = 'Two iterations of hemi as a smoke test'
    build_system = DummyBuildSystem()
    case = 'reframe/minihemi.case'

    @run_before('compile')
    def setup_case(self):
        self.executable = os.path.join(self.neko_build.install_dir, 'bin/neko')
