# Tridiagonal matrix solver

The **tridiagonal** class solves a tridiagonal matrix equation Ax = d, where A is an NxN matrix and x and d are vectors of length N. The matrix M is composed of three vectors of length N:

* a:  subdiagonal     (The first element is not used)
* b:  diagonal
* c:  superdiagonal   (The last element is not used)

## Constructor
*td* = **tridiagonal_t**(*a*, *b*, *c*) creates the initial matrix

## Methods
* **tridiatonal_t%solve**(*d*) Solves the matrix equation when given the right hand side d. The answer is returned in the input array.