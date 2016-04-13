function [V,R,e] = poly2basis(p,R)  
% function [V,R,e] = poly2basis(p,R)
%
% DESCRIPTION
%   Projects a vector of polynomials p onto the span of the monomials
%   contained in the vector R. 
%
% INPUTS 
%   p: 1-by-lp vector of polynomials.
%   R [Optional]: lr-by-1 basis of monomials. [ Default: R=monomials(p) ]
%
% OUTPUTS 
%   V: lr-by-lp matrix expressing the projection of the polynomial p 
%      on the monomials in R. The projection of p on to the span of
%      R is given by R'*V. 
%   R: Vector of monomials 
%   e: Difference between the input polynomial p and the projection
%      R'*V, i.e. e = p-R'*V.  If R contains all monomials in p then e=0.
%
% SYNTAX 
%   [V,R,e] = poly2basis(p,R);  
%
% EXAMPLE
%   pvar x1 x2;
%   p = [x1^2-9, 5*x1+3*x1*x2-4*x2^2];
%   [V,R,e] = poly2basis(p,monomials(p));
%   [V R]
%   p-R'*V
%
% See also monomials

% 1/28/08  PJS     Initial Coding

% Error checking
p = polynomial(p);
if nargin==1
    R = monomials(p);
end
if size(p,1)~=1
    error('p must be a 1-by-lp polynomial');
end
if ~ismonom(R) || size(R,2)~=1
    error('R must be a lr-by-1 vector of monomials');
end

% Get poly dimensions, number of terms, and number of variables
lp = length(p);
pdeg1 = p.degmat;
[ntp,nvp] = size(pdeg1);

lr = length(R); 
Rdeg1 = R.degmat;
ntr = size(Rdeg1,1);

% Create degree matrices of R/pi in a common set of variables
[allvars,tmp,idx] = unique( [p.varname; R.varname] ); 
nv = length(allvars);

pdeg = sparse(ntp,nv);
pdeg(:,idx(1:nvp)) = pdeg1;

Rdeg = sparse(ntr,nv);
[ridx,cidx] = find(R.coefficient);
Rdeg(:,idx(nvp+1:end)) = Rdeg1(ridx,:); 

% Write p(x) in terms of R
V = sparse(lr,lp); 
pcoef = p.coefficient;
[temp,ridx,pidx]=intersect(Rdeg,pdeg,'rows');   
for i1=1:lp    
   V(ridx,i1) = pcoef(pidx,i1); 
end

% V = sparse(lr,lp);
% pcoef = p.coef;
% for i1 = 1:lp
%     % Get i^th polynomial
%     picoef = pcoef(:,i1);
%     idx = find( picoef );
%     picoef = picoef(idx);
%     pideg = pdeg(idx,:);
%     
%     % Express pi in terms of monomials in R    
%     [temp,ridx,pidx]=intersect(Rdeg,pideg,'rows');   
%     V(ridx,i1) = picoef(pidx);    
% end

% Compute projection error
if nargout > 2
    e = p - R'*V;
end

