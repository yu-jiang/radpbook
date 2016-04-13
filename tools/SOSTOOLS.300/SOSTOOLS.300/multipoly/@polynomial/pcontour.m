function varargout = pcontour(p,v,domain,linespec,npts,var)
% function [C,h] = pcontour(p,v,domain,linespec,npts,var)
%
% DESCRIPTION 
%   Plots contours of p(x,y) at the contour values specified by the vector 
%   v. The contours are generated numerically by evaluating p on a grid of
%   values of x and y and then calling the CONTOUR function. 
%
% INPUTS 
%   p: 1-by-1 polynomial of two variables
%   v: N-by-1 vector of contour values (Default: v=1)
%   domain: 1-by-4 vector specifying the plotting domain,
%              [Xmin Xmax Ymin Ymax] 
%         (Default: domain = [-1 1 -1 1])
%   linespec: Color and linetype.  (Default:  linespec='b')
%   npts: 1-by-2 vector specifying the number of grid points along
%          each axis, [Num of X pts, Num of Y pts] 
%          (Default: npts = [100 100])
%   var: 1-by-2 vector of pvars specifying the x/y axis variables,
%              [Variable for X axis, Variable for Y axis]
%          (Default var = p.var)
%              
% OUTPUTS  
%   C,h: Contour matrix and contour handle object returned by CONTOUR
%  
% SYNTAX 
%  pcontour(p)
%  pcontour(p,v)
%  pcontour(p,v,domain)
%  pcontour(p,v,domain,linespec)
%  pcontour(p,v,domain,linespec,npts)
%  pcontour(p,v,domain,linespec,npts,var)
%  [C,h] = pcontour(p,v,domain,linespec,npts,var)
%
% EXAMPLE
%  pvar x y
%  p = (x-2)^2-(x-2)*y+y^2;
%  domain = [0 4 -2 2];
%  [C,h]=pcontour(p,[0.5 1 2],domain);
%  clabel(C,h);
%
% See also contour, clabel, pcontour3

% PJS 5/7/08   Initial coding
% PJS 11/23/10 Removed SOS call for domain

%-------------------------------------------------------------------
% Error Checking
%-------------------------------------------------------------------


%-------------------------------------------------------------------
% Set defaults 
%-------------------------------------------------------------------
if nargin==1
    v=[];
    domain=[];
    linespec=[];
    npts=[];
    var=[];
elseif nargin==2    
    domain=[];
    linespec=[];
    npts=[];
    var=[];
elseif nargin==3
    linespec=[];
    npts=[];
    var=[];
elseif nargin==4
    npts=[];
    var=[];
elseif nargin==5
    var=[];
end

% Default contour
if isempty(v)
    v=1;
end
lv = length(v);

% Default linespec
if isempty(linespec)
    linespec = 'b';
end

% Default npts
if isempty(npts)    
    Nx = 100;
    Ny = 100;
else
    Nx = npts(1);
    if length(npts)==1
        Ny = Nx;
    else
        Ny = npts(2);
    end
end

% Define variables as chars
if isempty(var)
    x = p.varname{1};
    y = p.varname{2};
elseif ispvar(var) || iscellstr(var)
    if ispvar(var)
        var = char(var);
    end
    x = var{1};
    y = var{2};
else
    error('var must be a vector of pvars');
end

% % Default domain
% if isempty(domain) && exist('sosopt.m')==2
%     % Compute default axis values using SOS optimization
%     xl = zeros(lv,1); xu = zeros(lv,1);
%     yl = zeros(lv,1); yu = zeros(lv,1);
%     for i1=1:lv
%         % Bounds
%         xl(i1) = -LOCALxub(p,-x,v(i1));
%         xu(i1) = LOCALxub(p,x,v(i1));
%         
%         yl(i1) = -LOCALxub(p,-y,v(i1));
%         yu(i1) = LOCALxub(p,y,v(i1));
%     end
%     xl = min(xl); xu = max(xu);
%     yl = min(yl); yu = max(yu);
%     
%     xpad = 0.1*(xu-xl);
%     domain = [0 0 0 0];
%     if xpad<=0
%         domain(1:2) = [xl-1 xl+1];
%     else
%         domain(1:2) = [xl-xpad, xu+xpad];
%     end
%     
%     ypad = 0.1*(yu-yl);
%     if ypad<=0
%         domain(3:4) = [yl-1 yl+1];
%     else
%         domain(3:4) = [yl-ypad, yu+ypad];
%     end
% elseif isempty(domain)
if isempty(domain)
    domain = [-1 1 -1 1];
end    
    
% Plot contour
xg = linspace(domain(1),domain(2),Nx);
yg = linspace(domain(3),domain(4),Ny);
[xg,yg] = meshgrid(xg,yg);
%pgrid = double(subs(p,{x,y},{xg,yg}));
pgrid = double(subs(p,{x; y},[xg(:)'; yg(:)']));
pgrid = reshape(pgrid,size(xg));
if lv==1
    % Single contour syntax for contour function
    v = [v v];
end

if nargout==0
    contour(xg,yg,pgrid,v,linespec);    
else
    [C,h]=contour(xg,yg,pgrid,v,linespec);     
    varargout = {C,h};
end
xlabel(x)
ylabel(y)
  


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
% [info,dopt] = sosopt(s,p.var,g);
% if isempty(dopt)
%     % Polynomial might be unbounded
%     xub = 1e3;
% else
%     xub = subs(g,dopt(:,1),dopt(:,2));
% end
