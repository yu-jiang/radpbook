function [A,B,f0] = plinearize(f,x,u,x0,u0)
% function [A,B,f0] = plinearize(f,x,u,x0,u0)
%
% DESCRIPTION 
%   This function linearizes the vector polynomial function f(x,u) about 
%   the trim point x=x0 and u=u0.  The linearizaztion is
%        f(x,u) = f(x0,u0) + A*z + B*w + Higher Order Terms
%   where z:=x-x0 and w:=u-u0 are the deviations from the trim values.
%
% INPUTS 
%   f: Vector field of polynomial system  (Ns-by-1 polynomial)
%   x: State  (Ns-by-1 vector of pvars)
%   u: Input  (Nu-by-1 vector of pvars)  [Optional]
%   x0: Trim state [Optional, Default: x0=0]
%   u0: Trim input [Optional, Default: u0=0]
%     
% OUTPUTS
%   A: State matrix
%   B: Input matrix
%   f0: f evaluated at (x0,u0)
%
% SYNTAX
%   [A,f0] = plinearize(f,x)
%   [A,f0] = plinearize(f,x,x0)
%   [A,B,f0] = plinearize(f,x,u)
%   [A,B,f0] = plinearize(f,x,u,x0)
%   [A,B,f0] = plinearize(f,x,u,x0,u0)
%
% EXAMPLE
%   pvar x1 x2 u;
%   x = [x1;x2];
%   f = [-2*x1+x2+x1^2-7; x1-3*x2+u+u^2+3];
%   x0 = [3;4];
%   u0 = 2;
%   [A,B,f0] = plinearize(f,x,u,x0,u0)
%
% See also jacobian, ptrim

% PJS 9/29/2009   Initial Coding
% PJS 10/05/2009  Update for linearizing f(x,u)


% Error Checking
if nargin==2
    %   [A,f0] = plinearize(f,x)
    u = [];
    x0 = [];
    u0 = [];
elseif nargin==3
    if ispvar(u)        
        %   [A,B,f0] = plinearize(f,x,u)
        x0 = [];
        u0 = [];        
    else
        %   [A,f0] = plinearize(f,x,x0)
        x0 = u;
        u = [];
        u0 = [];        
    end    
elseif nargin==4
    %   [A,B,f0] = plinearize(f,x,u,x0)
    u0 = [];
end

% Default trim values
if isempty(x0)
    Nx = length(x);
    x0 = zeros(Nx,1);
end
if isempty(u0)
    Nu = length(u);
    u0 = zeros(Nu,1);
end
   
% Evaluate function at trim
v = [x(:); u(:)];
v0 = [x0(:); u0(:)];

f0 = subs(f,v,v0);
f0 = double(f0);

% State matrix 
dfdx = jacobian(f,x);
A = subs(dfdx,v,v0);
A = double(A);

% Input matrix
if ~isempty(u)
    dfdu = jacobian(f,u);
    B = subs(dfdu,v,v0);
    B = double(B);
else
    B=f0;
end
