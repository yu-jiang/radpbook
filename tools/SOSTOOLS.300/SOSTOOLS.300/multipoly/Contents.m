% Multivariate Polynomial Toolbox
% Version 2.00, 23 November 2010.
%
% Creating polynomial objects
%    PVAR         - Construct a polynomial variable
%    MPVAR        - Construct a matrix or vector polynomial variable
%    POLYNOMIAL   - Construct a polynomial object
%    MONOMIALS    - Construct list of monomials
%    PDATAFIT     - Compute a polynomial least squares fit to data
%    PFUNCTIONFIT - Compute a polynomial least squares fit to a function
%
% Simulink:
%    POLYLIB.MDL  - Simulink block for polynomial objects
%
% Polynomial plotting:
%    PCONTOUR     - Plot 2d polynomial contours
%    PCONTOUR3    - Plot 3d polynomial contours
%
% Polynomial functions:
%    POLY2BASIS   - Project polynomial onto a basis of monomials
%    PLINEARIZE   - Linearize a vector polynomial function
%    PTRIM        - Find trim conditions for a polynomial dynamical system
%    PVOLUME      - Estimate the volume of a polynomial set
%    PSAMPLE      - Draw random samples from a polynomial set
%    PSIM         - Simulate a polynomial dynamical system
%    PPLANESIM    - Plot the phase plane for a polynomial dynamical system
%    INT          - Element-by-element integration of a polynomial
%    DIFF         - Element-by-element differentiation of a polynomial
%    JACOBIAN     - Compute Jacobian matrix of a polynomial vector
%    COLLECT      - Collect coefficients of specified variables
%    SUBS         - Symbolic substitution
%    CLEANPOLY    - Remove terms based on value of coefficient and degree
%
% Polynomial characteristics:
%    ISDOUBLE     - True for arrays of doubles
%    ISPVAR       - True for arrays of pvars
%    ISMONOM      - True for arrays of monomials
%    ISEMPTY      - True for empty monomials
%    ISEQUAL      - Element by element polynomial comparisons
%    SIZE         - Size of a polynomial matrix
%    LENGTH       - Length of a polynomial matrix
%    FIELDNAMES   - Get properties of a polynomial object
%
% Conversions:
%    P2S          - Convert from multipoly to symbolic toolbox
%    S2P          - Convert from symbolic toolbox to multipoly
%    DOUBLE       - Convert constant polynomial to a double
%    CHAR         - Converts a polynomial to its string representation.
%
% Overloaded arithmetic operations:
%    PLUS, +      - Add polynomials
%    MINUS, -     - Subtract polynomials
%    MTIMES, *    - Multiply polynomials
%    MPOWER, ^    - Power of a polynomial
%    HORZCAT, [,] - Horizontal concatentation of polynomials
%    VERTCAT, [;] - Vertical concatentation of polynomials
%    DIAG         - Diagonal poly matrices and diagonals of poly matrices
%    TRIL         - Extract lower triangular part of a polynomial matrix
%    TRIU         - Extract upper triangular part of a polynomial matrix
%    BLKDIAG      - Block diagonal concatenation of polynomial matrices
%    CTRANSPOSE, ' - Non-conjugate transpose of a polynomial
%    TRANSPOSE, .' - Non-conjugate transpose of a polynomial
%    RESHAPE      - Reshape a polynomial matrix
%    REPMAT       - Replicate and tile an array of polynomials
%    UPLUS        - Unary plus of a polynomial
%    UMINUS       - Unary minus of a polynomial
%    TIMES, .*    - Element-by-element multiply of polynomials
%    POWER, .^    - Element-by-element power of a polynomial
%    SUM          - Sum of the elements of a polynomial array
%    PROD         - Product of the elements of a polynomial array
%    TRACE        - Sum of the diagonal elements
%    DET          - Determinant of a polynomial matrix
%

