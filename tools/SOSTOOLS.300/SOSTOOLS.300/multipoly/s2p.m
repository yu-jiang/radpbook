function p=s2p(s)
% function p=s2p(s)
%
% DESCRIPTION
%   Converts from a symbolic toolbox polynomial to a multipoly polynomial.
%
% INPUTS
%   s: Polynomial created using the symbolic toolbox
%
% OUTPUTS
%   p: Polynomial created using the multipoly toolbox
%
% SYNTAX
%   p = s2p(s)
%
% See also p2s

% 1/30/2003: PJS  Initial Coding
% 11/23/2010 PJS  Updated for matrix polynomials and to use sym/eval

if ~isa(s,'sym')
    error('Input must be a symbolic toolbox object');
end

% Convert to variable precision arithmetic
s = vpa(s);

% Create pvars for each variable in s
vars = findsym(s);
if isempty(vars)
    nv = 0;
else
    vidx = [0 strfind(vars,',') length(vars)+1];
    nv = length(vidx)-1;
end
for i1 =1:nv
    varname = vars(vidx(i1)+1 : vidx(i1+1)-1);
    pvar(varname);
end

% Expand the polynomial 
se = expand(s);

% Evaluate in symbolic expression to convert to multipoly
p = eval(se);

