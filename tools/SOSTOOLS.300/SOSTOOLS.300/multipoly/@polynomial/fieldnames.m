function out=fieldnames(a)
% FIELDNAMES Get POLYNOMIAL object property names.
%  
% NAMES=FIELDNAMES(P) returns a cell array of strings containing the
% names of the properties associated with uncertain polynomial object, P.
%
% See also GET, SET, POLYNOMIAL
  
% 11/12/2002: PJS  Initial Coding
%    This function allows tab completion of a public data fields.
  
aprop = PVget(a,'PropNames');
out = aprop.GPropNames;
