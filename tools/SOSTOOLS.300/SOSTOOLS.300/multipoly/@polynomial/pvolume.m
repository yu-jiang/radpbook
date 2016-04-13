function [vol,volstd] = pvolume(p,v,domain,npts)
% function [vol,volstd] = pvolume(p,v,domain,npts)
%
% DESCRIPTION 
%   Estimate the volume contained in the set {x : p(x)<=v} using Monte
%   Carlo sampling.  npts are drawn uniformly from a hypercube and the
%   number of points, nin, contained in the set { x : p(x) <= v} is
%   counted.  The volume is estimated as vol = nin/npts.  An estimate
%   of the standard deviation of this volume is also computed.
%
% INPUTS 
%   p: 1-by-1 polynomial of n variables
%   v: scalar specifying the sublevel of the polynomial (Default: v=1)
%   domain: n-by-3 array specifying the sampling hybercube.  domain(i,1)
%           is a pvar in p and domain(i,2:3) specifies the min and max
%           values of the cube along the specified variable direction,
%              [X1, X1min, X1max; ...; Xn, Xnmin, Xnmax]
%         (Default: domain = [-1 1] along all variable directions)
%   npts: scalar specifying the number of sample points 
%          (Default: npts = 1e4)
%
% OUTPUTS  
%   vol: Volume estimate of { x : p(x)<= v}
%   stdvol: Standard deviation of the volume estimate.
%  
% SYNTAX 
%  pvolume(p)
%  pvolume(p,v)
%  pvolume(p,v,domain)
%  pvolume(p,v,domain,npts)
%  [vol,stdvol] = pvolume(p,v,domain,npts)
%
% EXAMPLE
%  pvar x1 x2
%  p = x1^2 + x2^2;
%  r = 2;
%  domain = [x1, -r, r; x2, -r, r];
%  [vol,stdvol] = pvolume(p,r^2,domain);
%  truevol = pi*r^2;
%  [truevol, vol]
%  [abs(truevol-vol) stdvol]

% PJS 5/13/08   Initial coding
% PJS 11/23/10  Updated std for volume est, removed SOS call for domain

%-------------------------------------------------------------------
% Error Checking
%-------------------------------------------------------------------

%-------------------------------------------------------------------
% Set defaults 
%-------------------------------------------------------------------
if nargin==1
    v=[];
    domain=[];
    npts=[];
elseif nargin==2    
    domain=[];
    npts=[];
elseif nargin==3
    npts=[];
end

% Default contour
if isempty(v)
    v=1;
end

% Default npts
if isempty(npts)    
    npts = 1e4;
end

% Default domain
% if isempty(domain) && exist('sosopt.m')==2    
%     var = polynomial(p.var);
%     nvar = length(var);
%     
%     xmin = zeros(1,nvar);
%     xmax = zeros(1,nvar);    
%     for i1=1:nvar
%         x = var(i1);
%         xmin(i1) = -LOCALxub(p,-x,v);
%         xmax(i1) = LOCALxub(p,x,v);                
%     end    
if isempty(domain)
    var = p.varname; 
    nvar = length(var);
    xmin = -ones(1,nvar);
    xmax = ones(1,nvar);
else
    nvar = size(domain,1);
    
    L.type = '()';
    L.subs = {':',1};
    var = subsref(domain,L);
    L.subs = {':',2};
    xmin = double( subsref(domain,L) )';
    L.subs = {':',3};
    xmax = double( subsref(domain,L) )';
end    

% Hypercube volume 
xdiff = xmax-xmin;
cubevol = prod(xdiff);

% Generate samples
xdiff = repmat(xdiff,[npts 1]);
xmin = repmat(xmin,[npts 1]);
xpts = xmin+xdiff.*rand(npts,nvar);
pv = double(subs(p,var,xpts'));

% Estimate Volume and Standard Deviation
% Ref: Monte Carlo Integration Wikipedia Entry
f = pv<=v;
fmean = sum(f)/npts;
fvar = sum( (f-fmean).^2 )/(npts-1);

vol = fmean*cubevol;
volvar = cubevol^2*fvar/npts;
volstd = sqrt(volvar);

% Old Code
% f = sum( pv <= v)/npts;
% vol = f*cubevol;
% volstd = sqrt( f*(1-f)/npts );

% %------------------------------------------------------------------------
% % Local function which computes an upper bound on:
% %       max x   s.t. p(x,y)=v
% % The bound is computed by relaxing the optimization to:
% %       min g    s.t. -x-h*(v-p)+g in SOS
% %------------------------------------------------------------------------
% function  xub = LOCALxub(p,x,v)
% 
% pvar g;
% h = polydecvar('c',1,'vec');
% 
% s=-x-h*(v-p)+g;
% [info,dopt,sossol] = sosopt(s,p.var,g);
% if ~info.feas
%     % Polynomial might be unbounded
%     xub = 1e3;
% else
%     xub = subs(g,dopt(:,1),dopt(:,2));
% end
