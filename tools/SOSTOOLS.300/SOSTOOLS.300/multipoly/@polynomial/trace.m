function b = trace(a)
% function B=trace(A)
%
% DESCRIPTION 
%   Sum of diagonal elements
%   
% INPUTS 
%   A: polynomial 
%
% OUTPUTS  
%   B: sum of the diagonal elements of A. 
%  
% SYNTAX 
%   B = trace(A);

% 6/7/2002: PJS  Initial Coding  
  
sza = size(a);
if sza(1)~=sza(2)
    error('Polynomial Matrix must be square');
elseif isempty(a)
    b = polynomial(0);
    return;
else
    % Grab diagonal elements
    acoef = a.coefficient; 
    idx = find(eye(sza));
    acoef = acoef(:,idx);

    % Sum diagonal elements
    coefficient = sum(acoef,2);
    b = a;
    b.coefficient = coefficient;
    b.matdim = [1 1];
    b = combine(b);
end
   