# Function library
The function library **functions** contains useful mathematical functions used frequently in many contexts

* *sp* = **softplus**(*s*,*x*) Softplus function sp(x) = ln(1 + exp(x-s)) that accepts scalar or vector values of *x*.
* *l* = **logistic**(*s*, *x*) Logistic function l(x) 1/(1 + exp(-(x-s)) that accepts scalar of vector values of *x*.
* *g* = **gaussian**(*x*, *x0*, *Delta*) Normalized Gaussain function g(x) = (Delta sqrt(2\*pi))^(-1) \* exp( (x-x0)^2/(2\*Delta^2) ) that accepts scalar or vector values of *x*.
* *r* = **rosenbrock**(*x*, *y*[, *a*, *b*]) Rosenbrock function r(x) = (a-x)^2 + b*(y-x^2)^2) that accepts scalar or vector values of *x* and *y*. Optional arguments are *a*=1 and *b*=100 by default.

