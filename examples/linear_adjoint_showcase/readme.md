# PowerIterations

We hav a split field, being driven by a baseflow and described by a field of
pertubations.

My goal ios to implement a simulation component which can compute the spectral
values of the flow. This way we should be able to describe the stability of the
system.
The spectral values are computed by the power iterations method.

## Power Iterations

The power iterations method is a method for computing the largest eigenvalue of
a matrix. The method is based on the fact that the largest eigenvalue of a
matrix is the limit of the ratio of the norm of the matrix to the norm of a
vector iteratively multiplied by the matrix.

For a given time step `t` we compute $\lambda_t$ based on the perturbation
field `u`:

$$
    \lambda_{t} = \frac{ \sum_{GLL} u_{t} \cdot \bar{u}_{t-1} }
        {\sum_{GLL} \bar{u}_{t-1} \cdot \bar{u}_{t-1}} \\
    \bar{u}_{t} = \frac{ u_{t} }{ ||u_{t}||_2 }\\
    ||u||_2 = \sqrt{ \frac{c}{V} \sum_{GLL} m (u_{t} \cdot u_{t}) }
$$

Where:
- $u_{t}$ is the perturbation field at time `t`.
- $\bar{u}_{t}$ is the normalized perturbation field at time `t`.
- $m$ is the mass matrix component for the gll points.
- $c$ is a constant which in Nek5000 is set to 0.5.
- $V$ is the volume of the domain.



## Implementation

This is already implemented for Nek5000 and we are going to use that as the
baseline for our work here.

- https://github.com/KTH-Nek5000/KTH_Framework
  - Examples/ext_cyl_ARN
  - Examples/ext_cyl_PWI
  - Toolbox/tools/tstpr/pwit/pwit.f
  - Toolbox/tools/tstpr/tstpr.f



# Random notes:
Compute the norm of the pertubation fields.

Nek5000 uses the following norm:
do il=1,ntotv
sum = sum + wt(il)*f1*( b1(il)*x1(il) + b2(il)*x2(il) + b3(il)*x3(il) )
 end do

Where:
      Harrison Nobis: Harrison Nobis said:
Harrison Nobis: a bit of Nek5000 notation:
if you want to take an integral of any scalar field in the domain, it's the pointwise product of the field and the mass matrix, then summed across the entire domain
so often people will use: glsc2(field1, bm1, n)
Harrison Nobis: in this case they're doing the sum locally on each processor, and then the last line glsum() does a global sum across all processors
Harrison Nobis: so by Neko notation, you would want to have:
uucoef%X_h%B + vvcoef%X_h%B + wwcoef%X_h%B
Harrison Nobis: I could be wrong with where B lives, but that's should be the mass matrix
Harrison Nobis: and then sum them everywhere




## Chat

### Harrison Nobis

ok here is the KTH toolbox https://github.com/KTH-Nek5000/KTH_Framework

The Arnoldi example I was working with is Examples/ext_cyl_ARN

The power iterations is a bit mislead, it's Examples/ext_cyl_PWI But in that
example we're doing optimal initial conditions, so we're direct-adjoint looping
we don't need to do that, we can just do power iterations on the direct, then
power iterations on the adjoint

Toolbox/tools/tstpr/pwit/pwit.f and Toolbox/tools/tstpr/tstpr.f

are the relevent modules for power iterations

### Tim Felle Olsen

Great, thanks Mate. I am going to look into that.

We should work on the linearized_NS branch right?

### Harrison Nobis

ok little snag, the baseflow is 2D and neko needs to read 3D files, so I'm going
to try convert it in matlab or so... but my optimism is draining rapidly!

### Tim Felle Olsen

### Harrison Nobis said:


ok little snag, the baseflow is 2D and neko needs to read 3D files, so I'm going
to try convert it in matlab or so... but my optimism is draining rapidly!

Oh crap..

### Harrison Nobis

few! Ok matlab came to the rescue, back on track

### Tim Felle Olsen

Hey Mate,

Two things.

Did you make a fork of the Neko repository and are you working on that. Or are
we working on the feature/linearized_NS branch? Did you test the
feature/linearized_NS branch currently pushed to Neko on a gpu? I tried to
execute the two examples (TS_channel/linear and TS_channel/nonlinear) but they
fail to run with an end time of 1s.
### Harrison Nobis

2) nothing on GPU's yet

1) I'm working locally on linearized_NS, and I'm almost ready to show you

Ok I finally got a case running... it's really awkward with the baseflow etc

but I did a super quick test and at least it can load the baseflow, start with
an initial random perturbation, and it seems to get advected by the baseflow

remember.... the whole point of this is to test the implementation of the
linearized NS convective term... so be SUUUUPPER critical of anything that looks
suspicious

But I guess you're having lunch now (I'm about to) so maybe we can talk after
lunch

the case is called DELETE_ext_cyl_PWI

give it a makeneko ext_cyl.f90

(also... we don't have a sponge at the end which is different to the Nek5000
example)

### Tim Felle Olsen

### Harrison Nobis said:


2) nothing on GPU's yet

Alright, could we signify that to the user with a better message than crashing
neko.

Maybe adding this in a strategically smart place?:

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error('The linearized solver does not support devices yet')
    end if
### Harrison Nobis said:


But I guess you're having lunch now (I'm about to) so maybe we can talk after
lunch

Sure lets talk a bit when you are back.

I can see that I am missing the fld file that we attempt to read in.

### Harrison Nobis

### Tim Felle Olsen said:


### Harrison Nobis said:


2) nothing on GPU's yet

Alright, could we signify that to the user with a better message than crashing
neko.

Maybe adding this in a strategically smart place?:

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error('The linearized solver does not support devices yet')
    end if

I love how forward thinking you are. Please make any changes and adjustments you
see fit

### Tim Felle Olsen said:


### Harrison Nobis said:


But I guess you're having lunch now (I'm about to) so maybe we can talk after
lunch

Sure lets talk a bit when you are back.

I can see that I am missing the fld file that we attempt to read in.

sorry!

harry.nek5000 harry0.f00000 harry0.f00001

I'm guessing github ignores the f files?

### Tim Felle Olsen

Yep, in general the fluid files are seen as outputs of Neko so they are ignored
to avoid bombarding the repository with massive files.

 12:53
### Harrison Nobis said:


harry.nek5000 harry0.f00000 harry0.f00001

Thanks Mate.

### Harrison Nobis said:


### Tim Felle Olsen said:


### Harrison Nobis said:


2) nothing on GPU's yet

Alright, could we signify that to the user with a better message than crashing
neko.

Maybe adding this in a strategically smart place?:

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call neko_error('The linearized solver does not support devices yet')
    end if

I love how forward thinking you are. Please make any changes and adjustments you
see fit

Well its not forward thinking as i just ran into the issue
:stuck_out_tongue_wink:

### Harrison Nobis

ok I'm done with lunch

I'm going to start implementing the adjoint operator now, do we want to have a
quick call first?

(Sorry I realized that "quick example case" took up the whole morning...
:grimacing: )

### Tim Felle Olsen

sure lets chat

### Harrison Nobis

im in your room

### Tim Felle Olsen

### Harrison Nobis said:


im in your room

that sound scary out of context :laughing:

### Harrison Nobis

eigenvalues_direct.txt eigenvalues_new_adjoint.txt eigenvalues_old_adjoint.txt

the first one comes from LNS The "new" adjoint has no integration by parts The
"old" adjoint has integration by parts + bdry

point is, they should all share the same spectrum

### Tim Felle Olsen

Great, thanks for the eigen value stuff.

What should be the output of the example you gave. I just get a field of all 0
entries when i load it up in paraview.

### Harrison Nobis

oh that's no good

I'm guessing since you're using GPU's you need to do a cheeky memcopy somewhere?

### Tim Felle Olsen

No no, I ran it with no GPU

### Harrison Nobis

did you compile with ext_cyl.f90 too? my first timestep is noise

image.png


### Tim Felle Olsen

Yep,. It will not run without it since you set the user_defined initial
condition.

### Harrison Nobis

and then it develops downstream

image.png


### Tim Felle Olsen

well at the end time, nothing is still there what kind of flowrates do you
expect to see after the 2s mark?

### Harrison Nobis

ooooooh at the end!

I didn't run it for long

but if it's decaying it means either I've fucked up the linear solver (which is
likely, but at least that's something to look at and debug) or I've fucked up
the baseflow (...also likely)

because in principle something should be growing

### Tim Felle Olsen

Well I get values for the flow magnitude in the range (0.0, 1e-8) which I assume
is just the precision we are working at.

### Harrison Nobis

but you're right, I can also see mine is decaying

the initial perturbation has amplitude 1.0e-08

Because it's assumed to be growing exponentially

but you're right that's well spotted..., something off

let me run for a few timesteps in nek5000 and see what happens

### Tim Felle Olsen

Glad im not crazy then.

I don't see the BC's anywhere. are those encoded in the mesh or how does that
work for this case?

Ah, I think i found it.

There is no inflow

its 0

### Harrison Nobis

they're supposed to be zero everywhere

### Tim Felle Olsen

Setting the

"inflow_condition": { "value": [1.0,0.0,0.0], } Fixed it

### Harrison Nobis

no it's still supposed to be zero everywhere...

### Tim Felle Olsen

How are we driving a flow with no inflow?

### Harrison Nobis

the baseflow drives the flow

### Tim Felle Olsen

I guess the baseflow is not applied then

### Harrison Nobis

OK I've just run it in Nek5000 and it's decaying aswell, so it seems ok

### Tim Felle Olsen

Alright.

### Harrison Nobis

hang on.... give me 1 sec to compare neko and nek5000

### Tim Felle Olsen

So the u,v,w that is being output from Neko is just the pertubation fields.

When we want to compute the Eigen values of the flow. is that just for the
pertubation field or is it on the full field (pertubation + base)

^Yes

