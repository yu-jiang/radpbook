function [g0,g,h] = collect(p,x)
% function [g0,g,h] = collect(p,x);
%
% DESCRIPTION
%   Collect p(x,y) into the form g0(x)+g(x)*h(y) where h(y) is a vector
%   of unique monomials in y.
%
% INPUTS
%   p: M-by-1 polynomial in variables x and y.
%   x: variables of p to collect into polynomials with coefficients given
%      by monomials in y. x can either be a polynomial vector or
%      a cell array of strings of variable names.
%
% OUTPUTS
%   g0: M-by-1 polynomial in x.
%   g: M-by-N vector of polynomials in x.
%   h: N-by-1 vector of monomials in y.
%
% SYNTAX
%   [g0,g,h] = collect(p,x);
%      g0, g, and h satisfy p(x,y) = g0(x)+ g(x)*h(y)
%
% EXAMPLE
%   pvar x1 x2 y1 y2;
%   p = 13+x1^2*y1-5*x1^2*y2^3+6*x1*x2*y1+8*x1;
%   x = [x1;x2];
%   [g0,g,h] = collect(p,x)
%   p-(g0+g*h)

% 1/29/08  PJS     Initial Coding
% 5/21/09  PJS     Updates for speed
% 10/27/10 PJS     Allow p to have dimension M-by-1

% XXX Error checking
szP = size(p);
if szP(2)>1
    error('Input polynomial p must be a column vector');
end
M = szP(1);

% Get list of "x" variables as a cell string
p = polynomial(p);
pcoef = p.coefficient;
if isa(x,'polynomial');
    x = x.varname;
end

% Rewrite degree matrix of p in common set of variables
[allvars,tmp,idx] = unique( [p.varname; x] );
nv = length(allvars);

[ntp,nvp] = size(p.degmat);
pdeg = sparse(ntp,nv);
pdeg(:,idx(1:nvp)) = p.degmat;

% Split out degree matrices associated with x and y variables
nx = length(x);
xidx = idx(end-nx+1:end);
xdegmat =  pdeg(:,xidx);
yidx = 1:nv;
yidx(xidx) = [];
y = allvars(yidx);

if isempty(yidx)
    % Column dim of h should equal szP(2), e.g. if p is 4-by-0
    g0 = p;
    g = polynomial(zeros(M,0));
    h = polynomial(zeros(0,szP(2))); 
    return;
else
    ydegmat =  pdeg(:,yidx);
end

% XXX PJS 12/2/09 The unique command below creates hdegmat that reverses
% the alphabetical order of y. Put y in reverse alphabetical order here so
% that the end result (after the unique command) is in alphabetical order.
y = flipud(y);
ydegmat = fliplr(ydegmat);

% Find unique list of monomials in y
[hdegmat,I,J]=unique(ydegmat,'rows');
lh = size(hdegmat,1);

% h, as constructed may have a constant monomial, hi(y)==1.
% Pull out the polynomial, g0(x), associated with this monomial.
g0idx = find( sum(hdegmat,2)== 0 );
chkval = 0; % skip validity check in polynomial constructor
if isempty(g0idx)
    g0 = polynomial(zeros(M,1));
    h = polynomial( speye(lh), hdegmat, y, [lh 1],chkval);
    
    gcoef = [];
    for i1=1:M
        tmp = sparse(1:length(J),J,pcoef(:,i1));
        gcoef = [gcoef tmp];
    end
    g = polynomial(gcoef,xdegmat,x,[lh M],chkval)';
    
    % --- OLD Code for 1-by-1 p
    % g = polynomial(sparse(1:length(J),J,p.coefficient),xdegmat,x,[M lh],chkval);
    % ---
    
    % Combine Terms
    g = combine(g);
else
    
    g0 = polynomial( pcoef(J==g0idx,:), xdegmat(J==g0idx,:),x,[M,1],chkval);
    
    idx = [1:(g0idx-1) (g0idx+1):lh];
    gcoef = [];
    for i1=1:M
        tmp = sparse(1:length(J),J,pcoef(:,i1));
        gcoef = [gcoef tmp(:,idx)];
    end
    g = polynomial(gcoef,xdegmat,x,[lh-1 M],chkval)';
    
    % --- OLD Code for 1-by-1 p
    %     idx = [1:(g0idx-1) (g0idx+1):lh];
    %     gcoef = sparse(1:length(J),J,p.coefficient);
    %     g0 = polynomial(gcoef(:,g0idx),xdegmat,x,[M 1],chkval);
    %     g = polynomial(gcoef(:,idx),xdegmat,x,[M lh-1],chkval);
    % ---
    
    hdegmat(g0idx,:) = [];
    h = polynomial( speye(lh-1), hdegmat, y, [lh-1 1],chkval);
    
    % Combine Terms
    g0 = combine(g0);
    g = combine(g);
end
