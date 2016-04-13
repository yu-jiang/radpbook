function [xin,xon]=psample(p,x,x0,N)
% function [xin,xon]=psample(p,x,x0,N)
%
% DESCRIPTION
%   This function draws random samples from a set described by a
%   single polynomial inequality:
%             S:={ x : p(x)<=0 }
%   A gas dynamics model is used to generate the random samples.  This
%   method requires an initial feasible point x0 in S.  The function also
%   assumes that S is closed and bounded.
%
% INPUTS
%   p: 1-by-1 polynomial of x used to describe the set S.
%   x: Nx-by-1 vector of pvars.  These are the variables in p.
%   x0: Initial point in the set S (Nx-by-1 double).  The values in x0
%      should correspond to the ordering of variables in x.
%   N: Number of random samples to generate. (default: N=1)
%
% OUTPUTS
%   xin: Nx-by-N matrix with each column specifying an element in S.
%   xon: Nx-by-N matrix with each column specifying an element on the
%       boundary of S, i.e. p(xon(:,i))==0 for each i.
%
% SYNTAX
%   [xin,xon]=psample(p,x,x0)
%   [xin,xon]=psample(p,x,x0,N)
%
% EXAMPLE
%   % Sample unit disk
%   pvar x1 x2;
%   x = [x1;x2];
%   p = x'*x-1;
%   [xin,xon]=psample(p,x,zeros(2,1),1e3);
%   plot(xon(1,:),xon(2,:),'rx'); hold on;
%   plot(xin(1,:),xin(2,:),'bo');hold off;
%   legend('Samples on Boundary','Samples in Interior')
%   axis equal;

% Coding
%   PJS     8/2/06  Initial Coding
%   Ufuk Topcu 11/16/2006 modified to keep the point on the boundary
%   PJS     1/23/08 First column is no longer x0
%   PJS     4/27/2009  Clean-up of polyset_sample

% Error checking
if nargin==3
    N=1;
end
if ~( isa(p,'polynomial') && all(size(p)==[1 1]) )
    error('p must be a 1-by-1 polynomial');
end
if ~( ispvar(x) && min(size(x))==1)
    error('x must be a vector of pvars');
end
if ~( min(size(x0))==1 && length(x)==length(x0) )
    error('Length of x must equal length of x0.');
end
if double(subs(p,x,x0)) > 0
    error('Initial point x0 is not in the polynomial set S.');
end

% Get polynomial data and reorder degmat to match the order of vars in x
Nx = length(x);
xv = char(x(:));
pdeg = p.degmat;
pcoef = p.coefficient;
Nt = size(pcoef,1);
pv = p.varname;
[tf,pidx]=ismember(xv,pv);
if length(pidx)~=length(pv) || any(tf==0)
    error('The variables in p must match those in x.');
end
pdeg = pdeg(:,pidx);
mxdeg = get(p,'maxdeg'); %mxdeg = p.maxdeg;

% Pre-compute binomial coefficients
b = cell(mxdeg+1,1);
b{1} = 1;
b{2} = [1 1];
for i1=3:(mxdeg+1)
    b{i1} = sum([b{i1-1} 0;0 b{i1-1}]);
end

% Compute samples via gas dynamics:
xin = zeros(Nx,N);
xon = zeros(Nx,N);
x = x0;

for i1=1:N
    % Generate random search direction
    v = randn(Nx,1);
    
    % Form p(x+t*v)  [ Code below is much faster than calling subs ]
    pt = zeros(Nt,mxdeg+1);
    for i2=1:Nt
        % Form i2^th term
        m = full(pcoef(i2));
        for i3=1:Nx
            % Form i3^th monomial in i2^th term, (x_i3+t*v_i3)^d
            d = pdeg(i2,i3);
            mi = b{d+1}.*(x(i3).^(0:d)).*(v(i3).^(d:-1:0));
            
            % Convolve i3^th monomial with current term
            m = conv(m,mi);
        end
        
        % Store i2^th term
        lm = length(m);
        pt(i2,end-(lm-1):end) =  m;
    end
    pt = sum(pt);
    
    % Compute L,U s.t. x+tv is in S:={x: p(x)>=0} for all t in [L,U].
    %pt(end) = pt(end);
    tcollide = roots(pt);
    ridx = find( abs(imag(tcollide)) <= 1e-8*abs(real(tcollide)) );
    tcollide = tcollide(ridx);
    
    [tmin,minidx]=min(abs(tcollide));
    if tmin <= 1e-8*max(abs(tcollide))
        tcollide(minidx) = [];
        if all(tcollide<0)
            % Point is on the boundary and x+t*v is outside for t>0
            U = 0;
            L = max( tcollide );
        elseif all(tcollide>0)
            % Point is on the boundary and x+t*v is outside for t<0
            U = min( tcollide );
            L = 0;
        else
            % Point is on the boundary but we need to determine
            % which direction is feasible
            U = min( tcollide( tcollide>0 ) );
            L = max( tcollide( tcollide<0 ) );
            
            if polyval(pt,L/2)<=0
                % t<0 is feasible
                U = 0;
            elseif polyval(pt,U/2)<=0
                % t>0 is feasible
                L = 0;
            else
                % Neither direction is feasible
                U = 0; L=0;
            end
        end
    else
        % Point is in the interior
        U = min( tcollide( tcollide>0 ) );
        if isempty(U);
            U = inf;
        end
        
        L = max( tcollide( tcollide<0 ) );
        if isempty(L)
            L = -inf;
        end
    end
    
    % Randomly choose t in [L,U] and jump to x+tv
    if isinf(U) && isinf(L)
        t = randn;
        t_on = sign(rand-1/2)*inf;
    elseif isinf(U)
        t = max(L,randn);
        t_on = L;
    elseif isinf(L)
        t = min(U,randn);
        t_on = U;
    else
        t = L+(U-L)*rand;
        if U==0
            t_on = L;
        elseif L==0
            t_on = U;
        else
            t_on = L+(U-L)*round(rand);
        end
    end
    
    xin(:,i1) = x+t*v;
    xon(:,i1) = x+t_on*v;
    x = xin(:,i1);
    %     if norm(xon(:,i1)-x0)<1e-4
    %         keyboard
    %     end
end



