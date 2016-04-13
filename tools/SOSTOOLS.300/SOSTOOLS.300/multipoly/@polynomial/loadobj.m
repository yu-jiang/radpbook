function b = loadobj(a)
% function b=loadobj(a)
%
% DESCRIPTION 
%   Load filter for polynomial objects. This function is called when
%   a polynomial object is loaded from a *.mat file. 
%
% INPUTS 
%   A: polynomial
%
% OUTPUTS  
%   B: polynomial
%  
% SYNTAX 
%   B = loadobj(A);

% 5/8/2008 PJS  Initial Coding to deal with update to coef matrix
% 11/16/2010 PJS Update for loading into new Matlab object format

% Get number of terms and size of polynomial
Nterms = size(a.degmat,1);
sza = a.matdim;

% Coef matrix used to be stored as a Nterms*sza(1)*sza(2) column 
% but is now stored as a Nterms-by-(sza(1)*sza(2)) matrix
b = a;
b.coefficient = reshape(b.coefficient,[Nterms sza(1)*sza(2)]);

% I ported multipoly to the new Matlab object format. This causes
% old multipoly objects to load as structures. The line calls the
% constructor so that they are returned as polynomials.
if ~isa(a,'polynomial')
   b = polynomial(b.coefficient, b.degmat, b.varname, b.matdim); 
end

