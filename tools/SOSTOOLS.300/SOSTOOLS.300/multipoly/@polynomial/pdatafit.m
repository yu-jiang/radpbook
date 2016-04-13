function [pfit,cfit,info] = pdatafit(p,x,Xdata,Ydata,W)
% function [pfit,cfit,info] = pdatafit(p,x,Xdata,Ydata,W)
%
% DESCRIPTION
%   This function finds the coefficients of a multivariate polynomial that
%   best fits given data in a least-squares cost.  The data is fit
%   with a linear combination of polynomial basis functions:
%      p(x,c) = c1*f1(x)+c2*f2(x) + ... + ck*fk(x)
%   where f1, f2, ..., fk are the polynomial basis functions. pdatafit
%   computes the coefficients c1, c2, ..., ck that minimize the fitting
%   error in a weighted squares cost:
%      min_c  sum_i (W(i)*e(i))^2
%   where e(i) is the fitting error of the i^th data point, i.e.
%   e(i) := p(Xdata(i,:),c) - Ydata(i).
%
% INPUTS
%   p: 1-by-1 polynomial.
%   x: Nx-by-1 vector of pvars that specifies the independent variables
%       in p.  All other variables in p are considered to be coefficients.
%   Xdata: Nx-by-Npts matrix of input data values.  The i^th row of Xdata
%       gives the data values associated with x(i).
%   Ydata: 1-by-Npts vector of output data values
%   W (Optional): 1-by-Npts weighting vector  [Default: W=ones(1,Npts)]
%
% OUTPUTS
%   pfit: Least-squares polynomial fit
%   cfit: Nc-by-2 cell array of the optimal coefficients.  The first
%         column contains the coefficients (as chars) and the second
%         column contains the optimal values.  The subs command can be
%         be used to replace the coefficients in any polynomial with
%         their optimal values, e.g. pfit = subs(p,cfit).
%   info: Data structure containing the matrices in the least squares
%         problem.  info has the fields A, b, cfit, W, e.  This gives
%         the data of the least squares problem in the form:
%                 min_c || diag(W)*(A*c-b) ||_2.
%         e = A*cfit-b is the fitting error.
%
% SYNTAX
%   [pfit,cfit,info] = pdatafit(p,x,Xdata,Ydata)
%   [pfit,cfit,info] = pdatafit(p,x,Xdata,Ydata,W)
%
% EXAMPLE
%  Xdata = linspace(100,200);
%  Ydata = 1./Xdata;
%  pvar c0 c1 c2 x;
%  p=c0+c1*x+c2*x^2;
%  [pfit,cfit,info] = pdatafit(p,x,Xdata,Ydata)
%  plot(Xdata,Ydata,'bx',Xdata,double(subs(pfit,x,Xdata)),'r--')
%  legend('1/X','pfit'); xlabel('x');
%
% See also pfunctionfit

% PJS 10/05/2009  Initial Coding

% Initialize weighting vector
Npts = length(Ydata);
if nargin==4
    W =ones(1,Npts);
end
W = W(:);

% Vectorize and get data dimensions
x = x(:);
Ydata = Ydata(:);
if any(size(Xdata)==[1 1])
    % Make sure Xdata lines up for single variable case
    % (If there are multiple variables then user must make sure
    %  dimensions of Xdata are as described in the help)
    Xdata = Xdata(:)';
end

% Collect p(x,c) into the form g0(x)+g1(x)*c1 + ... + gN(x)*cN
% where ci are the coefficients to be fit with least squares
[g0,g,c] = collect(p,x);
zerop = polynomial(0);
if isequal(c,zerop)
    pfit = p;
    cfit = cell(0);
    return;
elseif ~ispvar(c)
    error('Input polynomial p must be linear in the fitting coefficients');
end

% Solve least squares problem: min_c  || W*(p(Xdata,c)-Ydata) ||_2
Nc = length(c);
A = zeros(Npts,Nc);
L.type = '()';
for i1=1:Nc
    L.subs = {i1};
    gi = subsref(g,L);
    A(:,i1) = double(subs(gi,x,Xdata))';
end
b = Ydata - double(subs(g0,x,Xdata))';

Awt = repmat(W,[1 Nc]).*A;
bwt = W.*b;
cls = Awt\bwt;

% Create output variables including info data structure
%cfit = [char(c(:)) num2cell(cls)];
cfit = [c(:) cls];
pfit = subs(p,cfit);

info.A = A;
info.b = b;
info.cfit = cls;
info.W = W;
info.e = A*cls-b;
