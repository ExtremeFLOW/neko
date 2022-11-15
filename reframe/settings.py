site_configuration = {
    'systems': [
        {
            'name': 'dt',
            'descr': 'Dardel test nodes',
            'hostnames': ['dt\d.pdc.kth.se'],
            'modules_system': 'lmod',
            'partitions': [
                {
                    'name': 'cpu',
                    'descr': 'CPU on dt0',
                    'scheduler': 'squeue',
                    'launcher': 'srun',
                    'access': ['-w dt0', '--exclusive'],
                    'environs': ['PrgEnv-cray'],
                    'processor': {
                        'num_cpus': 256,
                        'num_cpus_per_core': 2,
                        'num_cpus_per_socket': 128,
                        'num_sockets': 2
                    }
                },
                {
                    'name': 'gpu',
                    'descr': 'GPUs on dt2',
                    'scheduler': 'squeue',
                    'launcher': 'srun',
                    'access': ['-w dt2', '--exclusive'],
                    'environs': ['PrgEnv-cray'],
                    'modules': ['rocm/rocm', 'craype-accel-amd-gfx908'],
                    'devices': [
                        {
                            'type': 'gpu',
                            'arch': 'amd',
                            'num_devices': 2
                        }
                    ],
                    'extras': {
                        'select_device': './rocm_select_gpu_device'
                    },
                    'variables': [
                        ['MPICH_GPU_SUPPORT_ENABLED', '1']
                    ]
                }

            ],
            'prefix': '/dt1/${USER}/reframe'
        },
        {
            'name': 'github-actions',
            'descr': 'Github Actions runner',
            'hostnames': ['*'],
            'partitions': [
                {
                    'name': 'cpu',
                    'scheduler': 'local',
                    'launcher': 'mpirun',
                    'max_jobs': 1,
                    'processor': {
                        'num_cpus': 2,
                        'num_cpus_per_core': 1,
                        'num_cpus_per_socket': 2,
                        'num_sockets': 1
                    },
                    'environs': ['default']
                }
            ]
        },
    ],
    'environments': [
        {
            'name': 'PrgEnv-cray',
            'modules': ['PrgEnv-cray'],
            'cc': 'cc',
            'cxx': 'CC',
            'ftn': 'ftn',
            'target_systems': ['dt']
        },
        {
            'name': 'default',
            'cc': 'mpicc',
            'cxx': 'mpicxx',
            'ftn': 'mpif90'
        },
    ],
    'logging': [
        {
            'handlers': [
                {
                    'type': 'stream',
                    'name': 'stdout',
                    'level': 'info',
                    'format': '%(message)s'
                },
                {
                    'type': 'file',
                    'level': 'debug',
                    'format': '[%(asctime)s] %(levelname)s: %(check_info)s: %(message)s',
                    'append': False
                }
            ],
            'handlers_perflog': [
                {
                    'type': 'filelog',
                    'prefix': '%(check_system)s/%(check_partition)s',
                    'level': 'info',
                    'format': (
                        '%(check_job_completion_time)s|reframe %(version)s|'
                        '%(check_info)s|jobid=%(check_jobid)s|'
                        '%(check_perf_var)s=%(check_perf_value)s|'
                        'ref=%(check_perf_ref)s '
                        '(l=%(check_perf_lower_thres)s, '
                        'u=%(check_perf_upper_thres)s)|'
                        '%(check_perf_unit)s'
                    ),
                    'append': True
                }
            ]
        }
    ],
    'general': [
        {
            'timestamp_dirs': '%FT%T',
        }
    ],
}
