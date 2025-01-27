# Filter {#filter}

\tableofcontents

## Filter base type

`filter_t` is a base type for filters on `field_t` by calling subroutine `filter_t%apply(F_out, F_in)`, where `F_out` and `F_in` are the output and the input fields, respectively.

The derived type of `filter_t` includes `PDE_filter_t`, `elementwise_filter_t`, ...

## PDE filter

The PDE filter is defined in the derived type `PDE_filter_t`. ...

## Elementwise filter

The elementwise filter is defined in the derived type `elementwise_filter_t`. Basing on the spectral element discretization, it is performed by mapping the field into hierarchical polynomials element by element (currently only Legendre-like polynomials are supported) and then attenuating the kernel for different orders. It could be constructed by calling subroutine `elementwise_filter_t%init(json_t, coef_t)`. 

After constructing the object, one should also set up `elementwise_filter_t%elementwise_filter_type` and `elementwise_filter_t%trnsfr`.

`elementwise_filter_type` is a string to specify the type of usage of the polynomial for the filter.

| Name                         | Polynomials |
| ---------------------------- | -------------------------- |
| `nonBoyd`                    | i-th order polynomials L_i(x) |
| `Boyd`                       | L_i(x) for i<2, L_i(x) - L_{i-2}(x) for i>=2|

The kernel is defined through an array `real(kind=rp), allocatable :: trnsfr(:)` whose size is the number of polynomials. The indices of `elementwise_filter_t%trnsfr` corresponds to the polynomial order in ascending order. `elementwise_filter_t%trnsfr` is initialized to be all `1` such that the filter has no effect by default unless further defined. Once `elementwise_filter_t%trnsfr` is defined, one should call the subroutine `elementwise_filter_t%build_1d()` to settle the settings. One could also refer `dynamic_smagorinsky.f90` for an example of using `elementwise_filter_t`.