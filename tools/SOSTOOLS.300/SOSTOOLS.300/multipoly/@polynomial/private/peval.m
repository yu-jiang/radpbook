function val = peval(x,coef,deg)
% function val = peval(x,coef,deg)
%
% DESCRIPTION
%   Evaluate a polynomial at x.
%   XXXX-Currently assumes that the polynomial is either a row or
%   column vector
%   XXXX-No error checking for speed
%
% INPUTS
%   x: nvars-by-npts matrix each column of which specifies a point
%        at which to evaluate the polynomial.
%   coef: nterms-by-lp coefficient matrix where lp is the length
%        of the polynomial
%   deg: nterms-by-nvars degree matrix
%
% OUTPUTS
%   val: lp-by-npts matrix where each column specifies the value of
%        the polynomial evaluated at the corresponding column of x.
%
% SYNTAX
%   val = peval(x,coef,deg)

% 6/15/06  PJS  Initial Coding

% Compute dimensions
[nterms,nvars]=size(deg);
npts = size(x,2);
lp = size(coef,2);
coef = full(coef);
deg = full(deg);

% Compute monomials
xs = shiftdim(x,-1);            % xs is 1-by-nvars-by-npts
xrep = repmat(xs,[nterms 1 1]); % xrep is nterms-by-nvars-by-npts
drep = repmat(deg,[1 1 npts]);  % drep is nterms-by-nvars-by-npts

% npow = xrep;
% idx = find(drep~=1);
% npow(idx) = xrep(idx).^drep(idx);

npow = xrep.^drep;
nmonom = prod(npow,2);          % nmonom is nterms-by-1-by-npts

% Evalute polynomial--mulitply by coefs and sum monomials
nmonomrep = repmat(nmonom,[1 lp 1]); % nmonomrep is nterms-by-lp-by-npts
coefrep = repmat(coef,[1 1 npts]);   % coefrep is nterms-by-lp-by-npts
val = sum(coefrep.*nmonomrep,1);     % val is 1-by-lp-by-npts

% Reshape val from 1-by-lp-by-npts to lp-by-npts
if npts==1
    val = val(:);
elseif lp==1
    val = val(:)';
else
    val = squeeze(val);
end

