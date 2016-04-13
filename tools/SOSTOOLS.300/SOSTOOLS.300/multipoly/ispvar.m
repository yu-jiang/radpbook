function b = ispvar(a)
% function B = ispvar(A)
%
% DESCRIPTION
%   Returns true for a polynomial variable or an array of polynomial
%     variables.
%
% INPUTS
%   A: polynomial
%
% OUTPUTS
%   B: true if A is a polynomial variable or an array of polynomial
%      variables.  Otherwise B is returned as false.
%
% SYNTAX
%   B = ispvar(A);

% 4/29/06    PJS  Initial Coding
% 1/29/2008: PJS  Updated documentation

if ~isa(a,'polynomial')
    b = false;
    return
end

acoef = a.coefficient;
adeg = a.degmat;
idx = find(acoef);

b = false;
if all(acoef(idx)==1) && all( sum( abs(acoef) ,1)==1) ...
        &&  all( sum(adeg,2) == 1)
    b = true;
end

% CODE TO SEARCH FOR PVARS ELEMENT BY ELEMENT
% 
% sza = size(a);
% if ~isa(a,'polynomial')
%     b = false(sza);
%     return
% end
% 
% acoef = a.coefficient;
% adeg = a.degmat;
% 
% % pvars are only one term
% idx1 = sum( acoef~=0 ,1)==1;
% 
% % pvars coef = 1;
% idx2 = sum( acoef ,1)==1;
% 
% % pvars degmat corresponds to a single var
% idx3 = find( idx1 & idx2 );
% [ridx,cidx] = find( acoef(:, idx3 ) );
% idx4 = sum(adeg(ridx,:),2) ~= 1;
% idx3( cidx(idx4) ) = [];
% 
% % Create logical output
% b = false(sza);
% b( idx3 ) = true;
