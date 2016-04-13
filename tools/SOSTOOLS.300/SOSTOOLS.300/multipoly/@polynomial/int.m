function b = int(a,x,L,U)
% function B = int(A,X,L,U)
%
% DESCRIPTION
%   Element-by-element integration of a polynomial with respect
%   to a single variable.
%
% INPUTS
%   A: polynomial
%   X: Scalar polynomial variable [Optional with default X = A.varname{1}]
%   L: Lower limit of definite integral
%   U: Upper limit of definite integral
%
% OUTPUTS
%   B: polynomial
%
% SYNTAX
%   B = int(A,X)
%     Indefinite integral of the polynomial, A, with respect to X.
%     X should be a polynomial variable or string.   Integration is done
%     element-by-element if A is a matrix.
%   B = int(A,X,L,U)
%     Definite integral of A with respect to X from lower limit L to
%     upper limit U.
%   B = int(A,X,[L U]);
%     Equivalent to B = diff(A,X,L,U)
%
% EXAMPLE
%   pvar x y z;
%   a = 2*x^3 - 2*x*z^2 + 5*y*z;
%   b = int(a,x)
%   diff(b,x)-a
%   c = int(a,[0 1])
%
% See also: diff, jacobian

% 12/6/2010 PJS  Initial Coding

% Process Inputs/ Error Checking
if nargin==1
    x = a.varname{1};
    L = []; U = [];
elseif nargin==2
    if isa(x,'double')
        % B = int(A,[L,U])
        L = x(1);
        U = x(2);
        x = a.varname{1};
    else
        % B = int(A,X)
        L = []; U = [];
    end
elseif nargin==3
    if isa(x,'double')
        % B = int(A,L,U)
        U = L;
        L = x;
        x = a.varname{1};
    else
        % B = int(A,X,[L U])
        U = L(2);
        L = L(1);
    end
elseif nargin~=4
    error(['Invalid syntax for the "int" command. ' ...
        'Type "help int" for more information.'])
end

if ispvar(x) && length(x)==1
    x = x.varname{1};
elseif ~ischar(x)
    error('X must be a single polynomial variable or a string');
end

% Get polynomial info about A
acoef = a.coefficient;
adeg = a.degmat;
avar = a.varname;
sza = a.matdim;
Naterms = size(acoef,1);

% Find variable we are differentiating with respect to.
varnumb = find( strcmp(avar,x) );

% Perform indefinite integral
szb = sza;
if isempty(varnumb)
    % int a(y) dx = a(y)*x
    bcoef = acoef;
    bdeg = [adeg ones(Naterms,1)];
    bvar = [avar; x];
else
    % int a(y)*x^n dx = a(y)*x^(n+1) / (n+1)
    bdeg = adeg;
    nplus1 = bdeg(:,varnumb)+1;
    bdeg(:,varnumb) = nplus1;
    
    bcoef = acoef;
    bcoef = lrscale( bcoef, 1./nplus1 , []);
    
    bvar = avar;
end
    
chkval = 0; % skip validity check
b = polynomial(bcoef,bdeg,bvar,szb,chkval);

% Evaluate definite integral
if ~isempty(L)
    bU = subs(b,{x},U);
    bL = subs(b,{x},L);
    b = bU-bL;
end

% Combine any common terms
b = combine(b);
