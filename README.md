# flib
A Fortran Library

## Standard library (stl.f90)
### String manipulation
* *lstr* = **lowercase**(*str*) Converts specified string *str* to lower case. Values inside quotation marks are ignored.
* *ustr* = **uppercase**(*str*) Converts specified string *str* to upper case. Values inside quotation marks are ignored.

### Searching
* *low* = **binary_search**(*arr*, *val*) Performs a guaranteed O(log N) binary search on a sorted array *arr*. Given the  provided value *val*, returns lower index *low* of the segment containting *val*. If *val=0*, the value is below the smallest value *arr(1)* in the array. If *val=N*, the value is above the largest value *arr(N)* in the array. Does not check that *arr* is sorted

### File manipulation
* *N* = **count_lines**(*fname*) Returns number of lines *N* in file *fname*.

### Interpolation
* *vq* = **linear_interp**(*x*, *v*, *xq*) Linear interpolation from a set of points *x* with values *v*. Returns value(s) *vq* at the query point(s) *xq*. Assumes *x* is monotonically increasing (does not check).
* *vq* = **bilinear_interp**(*x*, *y*, *v*, *xq*, *yq*) Bilinear interpolation from a set of points (*x*,*y*) with values v. Returns value(s) *vq* at the query point(s) (*xq*,*yq*). Assumes *x* and *y* are monotonically increasing (does not check).

### Wrappers
These functions convert from class(*) to intrinsic types. For example *i* = **integer_**(*cs*), where *i* is a scalar or array of integers and *cs* is a scalar or array of class(*)

* **integer_**
* **real_**
* **character_**
* **logical_**
* **complex_**

## Function library (functions.f90)
* *sp* = **softplus**(*s*,*x*) Softplus function sp(x) = ln(1 + exp(x-s)) that accepts scalar or vector values of *x*.
* *l* = **logistic**(*s*, *x*) Logistic function l(x) 1/(1 + exp(-(x-s)) that accepts scalar of vector values of *x*.
* *g* = **gaussian**(*x*, *x0*, *Delta*) Normalized Gaussain function g(x) = (Delta sqrt(2*pi))^-1 * exp( (x-x0)^2/(2*Delta^2) ) that accepts scalar or vector values of *x*.
