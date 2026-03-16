# Checkpoint File Concerns

## Summary

`chkp_file_t` currently has temporary read-state stored as type components:

- `chkp_Xh`
- `sim_Xh`
- `space_interp`
- `global_interp`
- `mesh2mesh`

That state appears to be used only during `chkp_file_read()` and `chkp_read_field()`.

## Main Concerns

1. Lifecycle ownership is unclear.

`chkp_file_read()` initializes interpolation helpers and also frees them again at
the end of the read. That makes the helpers effectively local scratch state,
even though they live on the file object.

2. The component-based design encourages destructor coupling.

Once `generic_file_t` gained a deferred `free()`, it became tempting to clean up
`space_interp` and `global_interp` in `chkp_file_t%free()`. That immediately ran
into lifecycle ambiguity, because the same helpers are already freed inside the
read path.

3. Lower-level `free()` routines are not obviously idempotent.

In particular:

- `src/sem/interpolation.f90`
- `src/global_interpolation/global_interpolation.f90`

These destructors free device resources guarded by `c_associated(...)`, but do
not consistently reset the C pointers afterward. A second `free()` call may
therefore attempt to release the same device allocation again.

4. The current shape of the API suggests the wrong ownership model.

If interpolation helpers are meant to be scratch state for one read, they should
probably be local variables in `chkp_file_read()` and passed explicitly into a
private helper routine, instead of living on `chkp_file_t`.

## What Was True Before

Before adding `generic_file_t%free()`, the successful checkpoint read path did
not obviously leak the interpolation helpers, because `chkp_file_read()` already
called:

- `this%global_interp%free()`
- `this%space_interp%free()`

at the end of the routine.

## Recommended Follow-up

1. Decide whether `chkp_file_t` should own persistent interpolation state at all.
2. If not, move interpolation objects to locals in `chkp_file_read()`.
3. If yes, make the lower-level `free()` routines idempotent and document the
   ownership/lifecycle clearly.
4. Only after that, consider making `chkp_file_t%free()` release those members.
