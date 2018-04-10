# Interpolation library
The interpolation library **interpolation** currently provides univariate interpolation. Each interpolation method can be performed using the classes, where the data *v(x)* are passed on construction. Interpolation is evaluated for real sample points *xq* passed as values or arrays and returned in the variable *vq* using the subroutine *interpolate*.
```
type(linear_interpolator_t) :: l
l = linear_interpolate_t(x, v)
call l%interpolate(xq, vq)
```
Convenience functions can also be used
```
vq = linear_interpolate(x, v, xq)
```

## Univariate interpolators
* **linear_interpolator_t**(*x*, *v*[, *low_bc*, *high_bc*])
* **cubic_spline_interpolator_t**(*x*, *v*[, *low_bc*, *high_bc*, *low_f*, *high_f*])
* **pchip_interpolator_t**(*x*, *v*)

## Convenience functions
* *vq* = **linear_interpolate**(*x*, *v*, *xq*[, *low_bc*, *high_bc*])
* *vq* = **cubic_spline_interpolate**(*x*, *v*, *xq*[, *low_bc*, *high_bc*, *low_f*, *high_f*])
* *vq* = **pchip_interpolate**(*x*, *v*, *xq*)