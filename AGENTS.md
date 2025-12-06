# Neko guide for AI agents

# What is Neko?
Neko is computational fluid dynamics (CFD) software based on the Spectral
Element Method. It primarily targets high-fidelity scale-resolving simulations
of turbulent flows.  

## Helping users setup new Neko simulation cases
Neko is a simulation software, and a common scenario is that the user will ask
you to help setup a simulation. Alternatively, they may want tips about specific
settings or and explanation about a setup they find in an existing case. The instructions under this header are specifically for these scenarios. We
will refer to a concrete simulation to as "the case". 

### General guidelines
Below are the main guidelines you should follow for retrieving information and
making decisions about how to setup a case.

- A simulation case typically lives in its own separate folder, which will call
  the case folder.
- In input to the simulation consists of at least two files
  - A file with the computational mesh. More instructions on meshing follow.
  - A JSON file, which is either a .json or a .case (most often the latter),
    containing the configuration of the various components of the simulation.
    This file is referred to as the **case file**.
- In many cases, an additional Fortran source file, referred to as the **user
  file** is part of the case setup. This allows advanced customization.
- The folder `examples` contains already setup cases to demonstrate Neko's
  functionality. The `examples/README.md` contains an overview of the cases in
  the different sub-folders, and each of them has a `README.md` of its own. If
  you are tasked wih creating a new case, it is a very good idea to base it on
  one of the examples instead of starting from scratch. Think of the examples as
  templates for your work. It is a very good idea to ask the user, whether there
  is a particular example, which would serve as a good starting point for the
  setup.
- When unsure, try to prompt the user for possible details. Most CFD simulations
  will be unique and require a lot of input to setup correctly. "Educated
  guesses" may often not work. If you make assumptions, you should explicitly
  tell the user about them.

### The case file
- Your definitive reference for the contents of the case file is the file
  `doc/pages/user-guide/case-file.md`. You should acknowledge that you got
  access to this file.
- Save the case files as the name of the case folder plus .json, and try to pretty-format as good as you can.

### The mesh file
- The general guide for meshing is `doc/pages/user-guide/meshing.md`.
- There are essentially four possible options.
  - The user provides you a mesh file as an .nmsh and you use that.
  - The user has a mesh in .re2 format and you can convert that to .nmsh by
    running `rea2nbin` on that file.
  - The user has a mesh in a different format and you should try to help them
    convert it first to .re2 and then to .nmsh. This is the most difficult
    scenario but you can try to provide guidance based on the meshing.md
    referenced above.
  - There is no mesh, but the geometry is a box. In this case you can generate
    the mesh yourself using the `genmeshbox` utility. The source code for this
    utility is in `contrib/genmeshbox/genmeshbox.f90`. Study it so that you have
    perfect understanding of how it works. Then, follow the user's instructions
    about the mesh: its element sizes (uniform or not), periodicity. Observe
    that for non-uniform element size distributions `genmeshbox` takes the
    location of element edges from files. You may therefore need to generate
    these files based on the user's specifications.
- **Always** run the utility `mesh_checker` on the generated or existing .nmsh and make sure it reports no errors.

### The user file
- The general guide for the user file is is `doc/pages/user-guide/user-file.md`.
- The file is just a Fortran module.
- Name the user file as the name of the case folder plus .f90.
- If it is possible to setup the case without a user file, you should always do
  that. Only create a user file if it is really necessary for the case setup.
- A **very good** idea is to look into the `examples/programming` directory,
  which contains examples for some basic things you can do in the user file.
  Look at the `README.md` in that subfolder for a discription of what is in the
  files, but all the .f90 there are heavily documented.
  Other examples may also contain user files.
- General development guidelines apply:
  - `doc/pages/developer-guide/code-style.md` for the Fortran code style.
  - `doc/pages/developer-guide/important_types.md` provides an overview of
    the most important types for Neko. It is a very good idea to look at their
    implementation in the respective .f90 files under th `src` structure.
- At the end always run `makeneko` on the generated user file and makes sure it
  compiles.
