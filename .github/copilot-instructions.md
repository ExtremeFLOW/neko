# Neko — code review instructions

This file provides LLM instructions for reviewing pull requests. 

**MOST IMPORTANT**: Follow the workflow steps and associated checklists below!

## General workflow

1. Read the AGENTS.md file in the root of the repository, it provides a
   high-level overview of the project.
2. Determine if you are being run in a GitHub environment, or locally. 
   If it is GitHub, you will 
   - Add your comments directly at relevant lines in the PR.
   - No need to write a summary.  
   
   If ran locally, you will make a neat listed summary of your comments, with 
   references to line numbers and files.
3. Conduct a review of Fortran code based on the criteria listed under "Fortran 
   code review" below. Stick to the criteria!
4. See if you find any *critical* issues outside of above criteria, including in
   device code. If you do, report them, but you must be very sure the issues are
   real.
4. Carefully read the `doc/AGENTS.md` file to understand rules related to 
   documentation and the associated code review checklist (section "Reviewing
   documentation"). Following that particular checklist is **mandatory**.
5. Conduct the review of the documentation based on the checklist.
 
# Fortran code review

Your review is limited to the criteria listed below. Do not give general
opinions about code style or design unless they violate the rules stated here.
You are more or less acting as an intelligent static analyser. A **very 
important** focus area of your analysis is memory management.

The list of criteria to check are:

- Any subroutine that creates a local `allocatable` must also manually
  `deallocate` it before returning.
- Any type that has `pointer` or `allocatable` components must have a routine
  that nullifies pointers and deallocates allocatables. In concrete types this
  routine must be called `free`; in abstract types with `deferred` procedures it
  may be called `free_base` or similar, always with a `free_` prefix.
- If a type has both a destructor (`free` or similar) and a constructor (`init`
  or similar), the constructor must begin by calling the destructor.
- Pay attention to object ownership when `pointer` or `allocatable` components
  are involved. Make sure ownership is correctly handled.
- Make sure types are fully instantiated in the `init` routines. If there are
  several such routines, makes sure they all do it. 