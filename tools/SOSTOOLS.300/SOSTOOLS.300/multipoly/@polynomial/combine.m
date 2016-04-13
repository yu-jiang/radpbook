function b = combine(a,chkval)
% function B = combine(A)
%
% DESCRIPTION 
%   Combine terms and eliminate unnecessary variables.
%   
% INPUTS 
%   A: polynomial 
%
% OUTPUTS  
%   B: polynomial, the result of reduction.
%  
% SYNTAX 
%   B=combine(A);
  
% 10/29/2002: PJS  Wrote this wrapper for PVcombine
% 5/20/2009:  PJS  Added chkval (=0 to skip validity check, else = 1)

if nargin==1
    % Default is to skip the validity check in polynomial
    chkval = 0;
end

% Promote a to polynomial 
%a = polynomial(a);

% Get polynomial fields/dimensions
adeg = a.degmat;
avar = a.varname;  
[nta,nva]=size(adeg);
sza = a.matdim;
nra = sza(1);
nca = sza(2);
acoef = a.coefficient; 

if isempty(a)  
    b = a;
    return;  
elseif nva==0
    % Add terms of a constant polynomial
    coefficient = sum(acoef,1);
    matdim = [nra nca];
    coefficient = reshape(coefficient, matdim);  
    b = polynomial(coefficient);
    return        
else  
    % Find unique set of variables
    if nva>1    
        [newdeg,uniquevar] = PVuniquevar(adeg,avar);
    else
        newdeg = adeg;
        uniquevar = avar;
    end

    % Eliminate zero coefs
    [ridx,cidx] = find( acoef );
    ridx = unique(ridx);
    acoef = acoef(ridx,:);
    newdeg = newdeg(ridx,:);
    
    % Find unique set of monomials
    if nta>1
        [newcoef,uniquedeg] = PVuniqueterm(acoef,newdeg,[nra nca]);
    else
        newcoef = acoef;
        uniquedeg = newdeg;
    end

    % Eliminate terms with zero coefficient uncovered by PVuniqueterm
    [ridx,cidx] = find( newcoef );
    ridx = unique(ridx);
    coefficient = newcoef(ridx,:);
    degmat = uniquedeg(ridx,:);
    if isempty(coefficient)
        coefficient = sparse(1,nra*nca);
        degmat = sparse(1,0);
    end

    % Eliminate variables that no longer exist
    [ridx,cidx] = find( degmat );
    cidx = unique(cidx);
    varname = uniquevar(cidx);
    degmat = degmat(:,cidx);
    nut = size(degmat,1);
    nuv = size(degmat,2);

    % Check if polynomial is constant
    if nuv==0
        coefficient = sum(coefficient,1);
        matdim = [nra nca];
        b = polynomial( reshape(coefficient,matdim) );
        return;
    end
    
    % Create the polynomial
    matdim = [nra nca];
    b = polynomial(coefficient,degmat,varname,matdim,chkval);
    return  
end







