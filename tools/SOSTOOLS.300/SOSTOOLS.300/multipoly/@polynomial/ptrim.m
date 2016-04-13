function [xt,ut,ft,ht,flg] = ptrim(f,x,u,x0,u0,h,opts)
% function [xt,ut,ft,ht,flg] = ptrim(f,x,u,x0,u0,h,opts)
%
% DESCRIPTION 
%   This function solves for trim states and inputs for the polynomial
%   dynamical system
%      dx/dt = f(x,u)
%   The trim values (xt,ut) satisfy f(xt,ut)=0.  FSOLVE is used to
%   solve these nonlinear equations.  Initial guesses for the trim
%   state/input can be passed to FSOLVE.  Additional equality 
%   constraints on the trim condition can be specified in the form 
%   h(x,u)=0 where h is a polynomial vector.
%
% INPUTS 
%   f: Vector field of polynomial system  (Nx-by-1 polynomial)
%   x: State  (Nx-by-1 vector of pvars)
%   u: Input  (Nu-by-1 vector of pvars)  
%   x0: Initial guess for trim state [Optional, Default: x0=0]
%   u0: Initial guess for trim input [Optional, Default: u0=0]
%   h: Equality constraints (Nh-by-1 polynomial) [Optional]
%   opts: Options for fsolve. See fsolve help [Optional]
%     
% OUTPUTS
%   xt: Trim state (Nx-by-1 vector)
%   ut: Trim input (Nu-by-1 vector)
%   ft: f evaluated at (xt,ut) (Nx-by-1 vector)
%   ht: h evaluated at (xt,ut) (Nh-by-1 vector)
%      If ptrim was successful finding a trim point then ft:=f(xt,ut)
%      and ht:=h(xt,ut) will both be equal to zero
%   flg: Exit flag returned by fsolve
%
% SYNTAX
%   [xt,ut,ft,ht,flg] = ptrim(f,x,u)
%   [xt,ut,ft,ht,flg] = ptrim(f,x,u,x0,u0)
%   [xt,ut,ft,ht,flg] = ptrim(f,x,u,x0,u0,h)
%   [xt,ut,ft,ht,flg] = ptrim(f,x,u,x0,u0,h,opts)
%
% EXAMPLE
%   pvar x1 x2 u;
%   x = [x1;x2];
%   f = [-2*x1+x2+x1^2-7; x1-3*x2+u+u^2+3];
%
%   % Find a trim condition
%   [xt,ut,ft] = ptrim(f,x,u)
%
%   % Find a trim condition with x2=4
%   h = x2-4;
%   x0 = []; u0 = [];
%   [xt,ut,ft,ht] = ptrim(f,x,u,x0,u0,h)
%
% See also fsolve, plinearize

% PJS 10/7/2009   Initial Coding

% Error Checking
if nargin==3
    x0 = [];
    u0 = [];
    h = [];
    opts = [];
elseif nargin ==4
    u0 = [];
    h = [];
    opts = [];        
elseif nargin==5
    h = [];
    opts = [];
elseif nargin==6
    opts = [];
end

% Set default initial starting guesses 
Nx = length(x);
if isempty(x0)
    x0 = zeros(Nx,1);
end
if isempty(u0)
    Nu = length(u);
    u0 = zeros(Nu,1);
end

% Set default fsolve options
if isempty(opts)
    opts = optimset('Display','Off');    
end

% Trim equations: f(xt,ut)=0 and h(xt,ut)=0
trimeq = [f;h];

% Solve for trim states/inputs, z:=[x;u]
fh = @(z) double(subs(trimeq,[x; u],z)); 
z0 = [x0(:); u0(:)];
[zt,fht,flg] = fsolve(fh,z0,opts);

xt = zt(1:Nx);
ut = zt(Nx+1:end);
ft = fht(1:Nx);
ht = fht(Nx+1:end);
