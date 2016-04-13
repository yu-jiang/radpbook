function b = PVset(a,property,value)
% function B = PVset(A,Property,Value)
%
% DESCRIPTION
%   Private set for polynomial objects.
%
% INPUTS
%   A: a polynomial object.
%   Property: polynomial property
%   Value: value of the property
%
% OUTPUTS
%   B: new polynomial object.
%
% SYNTAX
%   B = PVset(A,Property,Value);

% 6/7/2002: PJS  Initial Coding
% 1/31/03:  SP   Added 'matdim' as a setable object
% 3/30/03:  PJS  Removed 'matdim' as a setable object (use 'reshape')

switch property
    case 'coefficient'
        a.coefficient = value;
    case 'degmat'
        a.degmat = value;
    case 'varname'
        a.varname = value;
    otherwise
        error(['Property ' property ' is Invalid']);
end
b=a;






