# Neko Reframe test suite

Requires Reframe to run (`spack install reframe` or `pip install reframe`).

```
reframe -C settings.py -c checks.py -r --performance-report
```

Substitute `-r` for `-l` to see what tests would run.

The device backend is chosen for partitions with a device with `type: gpu`.
Otherwise the cpu backend is used.

## Options

Run only a specific test:

```
# 2 iterations of hemi
-n MiniHemi

# Small TGV (512 elements)
-n Tgv8

# Medium TGV (32768 elements)
-n Tgv32
```

Control which configurations are used with environment variables:

```
# Real precision
NEKO_REAL=sp,dp

# Fluid scheme
NEKO_SCHEME=plan4,pnpn
```

## How to add a new system

The `settings.py` file defines how to run the suite on different systems.
Reference: https://reframe-hpc.readthedocs.io/en/stable/config_reference.html.

CPU partitions require processor information which can be automatically
generated using `reframe --detect-host-topology`. GPU partitions require a
device with `type: gpu`.

Also update `valid_systems` and potentially `valid_prog_environs` in `checks.py`.
