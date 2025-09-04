# Extending neko {#extending}

In addition to compiling the user file, `makeneko` can also compile extra `.f90`
files containing Fortran modules and `.cu/.hip` files containing CUDA/HIP
kernels. A common use of this is to separate out functionality needed by the
user file into individual source files. However, you can also use this feature
to extend Neko's core functionality.

Specifically, you can write modules that implement custom components such as LES
models, wall models, or simulation components. These modules can then be
configured directly via the case file's JSON input—just like Neko's built-in
components.

Of course, you could also achieve this by modifying Neko's source code directly.
But that comes with some drawbacks:

- You need to manually integrate your new files into Neko's build system.
- Updating Neko via Git becomes more complex, as you'll need to resolve merge
  conflicts between your changes and upstream updates.
- To avoid that hassle, you would have to open a pull request and merge your
  changes into Neko—though this may not be your intention. Or it might be to
  small in scope to qualify for being included.



Extending Neko with `makeneko` is quite easy, but not everything in the code
supports this. Here is a list of things you can implement.

- LES models (`les_model_t` descendants).
- Wall models (`wall_model_t` descendants).
- Simulation components (`simulation_component_t` descendants).
- Source terms (`source_term_t` descendants).
- Point zones (`point_zone_t` descendants).

This list will hopefully be extended later. Notable omissions are scalar and
fluid schemes, so currently you cannot add new solvers like this. 

To implement a new type, the easiest thing is to start by copying over the
`.f90` of an already existing type, renaming things inside and then adding the
functionality you want (e.g. the particular form of the source term). When you
are done, two special subroutines have to be added to you module. To make
explaining easier, let us assume that your module is called `mymodule` and the
new type `mytype`.

- All the types in the list above have a routine called `register_*`, where the
 `*` is replaced by the name of the type, for example, `register_source_term`.
- The also provide a procedure interface called `*_allocate`, e.g. 
  `source_term_allocate`.

Both of these routines should be brought into the the new modules with `use`
statements. Then, two routines need to be defined in the module.

- One is just an allocator for our new type. The name of the routine is 
  arbitrary.

  ```fortran
  subroutine mytype_allocate(obj)
    class(source_term_t), allocatable, intent(inout) :: obj

    ! The only thing the routine does is allocate a polymorphic object to our
    ! type.
    allocate(mytype_t::obj)
  end subroutine 
  ```
- The second routine uses the one above to register the type with Neko. It *has*
  to be called `mymodule_register_types`, where `mymodule` coincides with the 
  name of our module.

  ```fortran
  subroutine mymodel_register_types()
    ! Just a helper variable, not that we use the *_allocate routine, here
    ! for the source_term_t type, as an example
    procedure(source_term_allocate), pointer :: allocator

    ! Point to our allocator
    allocator => my_source_term_allocate

    ! Based on this the name of the source term will be "mytype",
    ! This is what you set in the JSON file, in the type keyword
    call register_source_term("mytype", allocator)
  end subroutine 
  ```

  Note that the `*_register_types` routine can register maybe types, not
  necessarily just one.

For custom device kernels, `mymodule` must define a C interface to a CUDA/HIP
routine that launches the kernel.

```fortran
interface
  subroutine device_kernel(a_d, n) &
        bind(c, name = 'device_kernel')
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr
    type(c_ptr), value :: a_d
    integer(c_int) :: n
  end subroutine device_kernel
end interface
```

Furthermore, the CUDA/HIP file must allow for C linkage, hence the routine
`device_kernel` must be inside an `extern "C"` block.

```cpp
extern "C" {
  void device_kernel(void *a, int *n) {
    /* Launch device kernel here */
  }
}
```

After compiling with `makeneko`, you can select your type in the JSON in the 
appropriate place.

One might wonder whether it's necessary to use this functionality—for example, 
for source terms—when it's possible to simply use a corresponding user routine 
in the user file. However, implementing a proper custom type has several 
advantages:

- You can potentially submit it to Neko via a pull request, if that's your goal.
- It's easier to distribute. For example, you can include the type in your own 
  Git repository, making it readily accessible to collaborators. While you could 
  also share a user file, combining multiple user files is error-prone, as they 
  often rely on module-level variables. Consider the difficulty of stitching 
  together two or three custom source terms in a single user routine without 
  causing conflicts.
- You can configure the type directly in the JSON input. This is especially 
  useful for parameter studies, where you may want to adjust values without 
  modifying code or recompiling.

That said, writing a custom type does take more effort than simply filling in 
the user routine.

