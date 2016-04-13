function [flag,errormsg] = PVisvalid(a)
% function [Flag,Errormsg] = PVisvalid(a)
%
% DESCRIPTION
%   (INTERNAL FUNCTION)
%   Check if a is a valid structure for a polynomial object.
%
% INPUTS
%   a: polynomial object. The validity of this structure will be checked.
%
% OUTPUTS
%   Flag: Flag=1 if structure is valid and 0 otherwise.
%   ErrorMsg: if Flag=0, this is a message indicating the
%      location and type of error.
%
% SYNTAX
%   [Flag,Errormsg] = PVisvalid(a);

% 6/9/2002: PJS  Initial Coding
% 10/22/2010: PJS  PVisvalid takes poly object in new TMW OO format

% Initialize Outputs
flag = 1;
errormsg = [];

% Check for errors in the structure
% 10/22/2010: In TMW new OO format, a is already a poly object
if ~isa(a,'polynomial')
    flag = 0;
    errormsg = 'Input must be a structure';
    return;
end

% Grab fields
adim = a.matdim;
acoef = a.coefficient;
adeg = a.degmat;
avar = a.varname;

% Empty polynomials are valid
if isempty(a.coefficient) && isempty(a.degmat) && ...
        isempty(a.varname) && all(adim==[0 0])
    flag = 1;
    errormsg = [];
    return;
end

% Check matdim: Should be 1x2 vector of integers
if ~( isa(adim,'double') && isreal(adim) && ndims(adim)==2 && ...
        all(floor(adim(:))==ceil(adim(:))) && all(size(adim)==[1 2]) ...
        && all(adim(:) >= 0) )
    flag = 0;
    errormsg='Matdim should be a 1x2 matrix of natural numbers';
    return;
end

% Check coefficient: Should be a real array.
if ~( isa(adim,'double') && isreal(acoef) && ndims(acoef)==2 )
    flag = 0;
    errormsg = 'Coefficient must be a real vector.';
    return;
end

% Check degmat:  Should be a matrix of natural numbers.
adegnz = nonzeros(adeg);
if ~( isa(adeg,'double') && isreal(adeg) && ndims(adeg)==2 && ...
        all(floor(adegnz)==ceil(adegnz)) && all(adegnz >= 0) )
    flag = 0;
    errormsg='Degmat should be a matrix of nautral numbers';
    return;
end

% Check size of Coefficient
nr = adim(1);
nc = adim(2);
nt = size(adeg,1);
if ~(  size(acoef,1)==nt && size(acoef,2)==(nr*nc) )
    flag = 0;
    errormsg='Size of Coefficient must be T x (N*M).';
    return;
end

% Check varname: Should be column cell arrays of strings with
% length == nv.
nv = size(adeg,2);
if length(avar) ~= nv;
    flag = 0;
    errormsg=['Number of variables implied by Degmat and Varname are ' ...
        'inconsistent.'];
    return;
end
if (nv>0) && ~( ndims(avar)==2 && iscellstr(avar) && size(avar,2) == 1 )
    flag = 0;
    errormsg='Varname should be a column cell array of strings';
    return;
end
if (nv==0) && ~( isempty(avar) && iscell(avar) )
    flag = 0;
    errormsg='Varname should be a column cell array of strings';
    return;
end
