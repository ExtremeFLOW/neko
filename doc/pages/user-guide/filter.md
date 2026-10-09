# Filtering {#filter}

\tableofcontents

## Filter base type {#filter_base_type}

`filter_t` is a base type for filters on `field_t` by calling subroutine `filter_t%apply(F_out, F_in)`, where `F_out` and `F_in` are the output and the input fields, respectively.

The derived type of `filter_t` includes `PDE_filter_t`, `elementwise_filter_t`, ...

## PDE filter {#filter_pde}

The PDE filter is defined in the derived type `PDE_filter_t`.
The filter is `based on the work of   [B. S. Lazarov, O. Sigmund]( https://doi.org/10.1002/nme.3072)
that solves a Helmholtz-type differential equation to provide smoothing. The
equation has the form
\f[
    -r^2 \nabla^2 X_\text{out} + X_\text{out} = X_\text{in},
\f]
subject to Neumann boundary conditions.
The filter has the
following input parameters:


| Name | Description  | Admissible values | Default value |
|------|--------------|-------------------|---------------|
| `r` | \f$r\f$ is the above equation. | Real | - |
| `tol`| The desired tolerance used when solving the system. | Real | `0.0000000001` |
| `max_iter` | Maximum number of iterations when solving the system. | Integer | `200` |
| `solver` | Numerical solver used to solve the system. | `cg`,`gmres`, `gmres` | `cg` |
| `preconditioner` | Pre-conditioner used to solve the system. | `ident`, `hsmg`, `jacobi` | `jacobi`  |



## Elementwise filter {#filter_elementwise}

The elementwise filter is defined in the derived type `elementwise_filter_t`. Basing on the spectral element discretization, it is performed by mapping the field into hierarchical polynomials element by element (currently only Legendre-like polynomials are supported) and then attenuating the kernel for different orders. One could directly set it up through the json file if the filter is needed in any section.

The user can also modify `elementwise_filter_type` and `transfer_function` according to their need.

`elementwise_filter_type` is an optional string to specify the type of usage of the polynomial for the filter.

| Name                         | Polynomials |
| ---------------------------- | -------------------------- |
| `nonBoyd`                    | i-th order polynomials L_i(x) |
| `Boyd`                       | L_i(x) for i<2, L_i(x) - L_{i-2}(x) for i>=2|

The kernel is defined through an optional array `transfer_function` whose size is the number of polynomials. The indices of `transfer_function` corresponds to the polynomial order in ascending order. `transfer_function` is initialized to be all `1` such that the filter has no effect by default unless further defined.

One could set up the elementwise_filter in the following way for polynomial order `7`.
~~~~~~~~~~~~~~~{.json}
"filter": {
    "type": "elementwise",
    "elementwise_filter_type": "Boyd",
    "transfer_function": [1,1,1,0.7,0.3,0,0,0]
}
~~~~~~~~~~~~~~~
One could also refer to `dynamic_smagorinsky.f90` for an example of using `elementwise_filter_t` in the neko code.

A list given in `transfer_function` is tied to one polynomial order and must
be changed with it. As an alternative, the optional `transfer` object builds
the transfer function from parameters given on the normalised mode number
\f$ k/N \f$, so the same case file works at every order. With \f$ N \f$ the
polynomial order, the modes up to \f$ k_c = \max(2, \mathrm{nint}(c\, N)) \f$
are kept and the higher ones are damped as
\f[
   \sigma(k) = \exp\left(-\alpha \left(\frac{k - k_c}{N - k_c}\right)^p\right),
   \qquad k > k_c ,
\f]
where \f$ c \f$ is `cutoff`, \f$ \alpha \f$ is `strength` and \f$ p \f$ is
`order`. The modes 0, 1 and 2 are never damped. In the Boyd basis they are the
only modes with a nonzero element integral, so the filter preserves the
element means of the filtered field. Give either `transfer_function` or
`transfer`, not both.

| Name                | Description                                                      | Admissible values | Default value |
| ------------------- | ---------------------------------------------------------------- | ----------------- | ------------- |
| `transfer.cutoff`   | Fraction of the spectrum that is kept.                           | Real in [0, 1]    | `0.4`         |
| `transfer.strength` | Damping exponent \f$ \alpha \f$ at the highest mode.             | Real, \f$ \ge 0 \f$ | `3.5`       |
| `transfer.order`    | Shape of the roll-off, 1 for exponential, 2 for Gaussian.        | Real, \f$ > 0 \f$ | `1.0`         |

The defaults give the following transfer functions.

| Order | Transfer function |
| ----- | ----------------- |
| 5     | 1, 1, 1, 0.31, 0.10, 0.03 |
| 7     | 1, 1, 1, 1, 0.42, 0.17, 0.07, 0.03 |
| 9     | 1, 1, 1, 1, 1, 0.50, 0.25, 0.12, 0.06, 0.03 |

~~~~~~~~~~~~~~~{.json}
"filter": {
    "type": "elementwise",
    "elementwise_filter_type": "Boyd",
    "transfer": {
        "cutoff": 0.4,
        "strength": 3.5,
        "order": 1
    }
}
~~~~~~~~~~~~~~~

## High-pass filter relaxation source term {#filter_hpfrt}

The high-pass filter relaxation source term uses the elementwise filter as a
low-pass filter and adds the high-pass residual to the right hand side. It can
be used as a source term for the fluid momentum equation or for a scalar
transport equation. See [Source terms](@ref case-file_fluid-source-term) for
the case-file setup.

Let \f$ q \f$ denote one velocity component or a scalar field. The source term
added to the equation is
\f[
   f_\text{hpfrt}(q) = -\chi (I - F) q,
\f]
where \f$ F \f$ is the tensor-product elementwise low-pass filter and \f$ \chi \f$ is the filter weight.
The negative sign makes the term dissipative for the high-pass part of the
field.

The one-dimensional low-pass filter is defined by a modal transfer function.
For `filter_modes` \f$ = k_f \f$ and \f$ n = lx \f$, let
\f[
   k_0 = n - k_f.
\f]
The transfer function used by the implementation is
\f[
\sigma_k =
\begin{cases}
1, & k \leq k_0,\\
1 - \left(\frac{k-k_0}{k_f}\right)^2, & k_0 < k \leq n.
\end{cases}
\f]
Thus the highest `filter_modes` modes are damped with a quadratic ramp, and the
highest mode is removed from the low-pass-filtered signal. The tensor-product
filter \f$ F \f$ is then applied element by element in each coordinate
direction.

For the momentum equation, the source term is applied componentwise to
\f$ u \f$, \f$ v \f$, and \f$ w \f$. For a scalar equation, it is applied to
the scalar field owned by that scalar scheme.
