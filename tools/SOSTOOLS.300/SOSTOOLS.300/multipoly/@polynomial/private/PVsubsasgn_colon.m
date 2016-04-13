function b = PVsubsasgn_colon(a,L,RHS)
% function B =  PVsubsasgn_colon(A,L,RHS)
%
% DESCRIPTION
%   (INTERNAL FUNCTION)
%   (:)-subsassign for polynomial objects.
%
% INPUTS
%   A: polynomial
%   L: a structure array with the fields:
%    type -- string containing '()'
%    subs -- Cell array or string containing ':'
%   RHS: Value to be assigned.
%
% OUTPUTS
%   B: object after subsassignment
%
% SYNTAX
%   B =  subsasgn(A,L,RHS)

% 6/9/2002: PJS  Initial Coding


% Get info about polynomials
b = a;
temp = RHS;
[nrb,ncb]=size(b);
[nrt,nct]=size(temp);

if nrt*nct == nrb*ncb
    if nrt*nct==0
        b = a;
    else
        b = reshape(temp,size(b));
    end
elseif nrt==0 || nct==0
    b = polynomial;
elseif all([nrt nct]==[1 1])
    if isempty(b)
        b = a;
    else
        b = temp;
        bcoef = b.coefficient;
        b.coefficient = repmat(bcoef,1,nrb*ncb);
        b.matdim = [nrb ncb];
    end
else
    error(['In an assignment  A(:) = B, the number of' ...
        ' elements in A and B must be the same.']);
end

