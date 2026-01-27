# Filtering {#filter}

\tableofcontents

## Filter base type

`filter_t` is a base type for filters on `field_t` by calling subroutine `filter_t%apply(F_out, F_in)`, where `F_out` and `F_in` are the output and the input fields, respectively.

The derived type of `filter_t` includes `PDE_filter_t`, `elementwise_filter_t`, ...

## PDE filter

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



## Elementwise filter

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
