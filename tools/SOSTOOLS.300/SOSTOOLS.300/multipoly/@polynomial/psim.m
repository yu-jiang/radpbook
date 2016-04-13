function [xtraj,xconv]=psim(f,x,x0mat,tfinal,event_params,opts)
% function [xtraj,xconv]=psim(f,x,x0,tfinal,event_params,opts)
%
% DESCRIPTION
%   Simulates non-autonomous polynomial systems of the form:
%               dx/dt = f(x),   x(t) = x0
%
% INPUTS
%   f: Vector field  (Ns-by-1 polynomial)
%   x: State  (Ns-by-1 vector of pvars)
%   x0: Initial Conditions (Ns-by-N0 array of doubles)
%   tfinal: Final simulation time unless the simulation is terminated
%     by one of the event parameters  (scalar)
%   event_params (Optional):  Event parameters for stopping the
%     simulation.  This is a structure with the following fields:
%      *nbig: Terminate if norm(x) is greater than nbig*norm(x0)
%      *nsmall: Terminate if norm(x) is less than nsmall*norm(x0)
%      *xbig: Terminate if any abs(x(i)) is greater than xbig(i)
%      *xsmall: Terminate if all abs(x(i)) are less than xsmall(i)
%      *funchandle:  Handle to a user specified event function.
%      *Additional fields can be used to pass parameter data to the
%          user defined event function
%      (Default: nbig = 1e6, nsmall = 1e-6, xbig =0, xsmall=0,
%          funchandle = [])
%   opts (Optional): Options structure passed to ODE solver. See
%       odeset and odeget for more details.  opts can have the additonal
%       field 'Solver' to specify the ode solver.  The 'Solver' field
%       can be ode45 or ode15s.
%
% OUTPUTS
%   xtraj: N0-by-2 cell array with the i^th row containing the simulation
%      results starting from x0(:,i).  xtraj{i,1} is an Nt-by-1 vector
%      of simulation times and xtraj{i,2} is an Nt-by-Ns matrix of
%      the state trajectories.
%   xconv: N0-by-1 logical vector with the i^th element = 1 if the
%      corresponding trajectory converged to the origin and = 0 otherwise.
%      A trajectory is considered to have converged to the origin if
%      either the nsmall or xsmall event occured.
%
% SYNTAX
%   [xtraj,xconv]=psim(f,x,x0,tfinal)
%   [xtraj,xconv]=psim(f,x,x0,tfinal,event_params)
%   [xtraj,xconv]=psim(f,x,x0,tfinal,event_params,opts)
%
% See also:
%  ODE solvers: ode45, ode23, ode113, ode15s, ode23s, ode23t, ode23tb
%  Options handling: odeset, odeget

% PJS 6/7/06   Initial coding of sim_rsys for rational systems
% PJS 4/22/09  Simplified sim_rsys for simulation of polynomial systems


% XXX Accept function handles in place of polys

% Error Checking
error(nargchk(4, 6, nargin, 'struct'))
if nargin==4
    event_params = [];
    opts = [];
elseif nargin==5
    opts = [];
end

if ( isa(f,'polynomial') && min(size(f))==1 )
    Ns = length(f);
else
    error('f must be a polynomial vector.');
end
if ~( ispvar(x) && min(size(x))==1 && length(x)==Ns )
    error('x must be a vector of pvars of the same length as f.');
end
if ( isa(x0mat,'double') && ndims(x0mat)==2 && size(x0mat,1)==Ns )
    N0 = size(x0mat,2);
else
    error('x0 must be a double array of dimensions Ns-by-N0');
end
if ~( isa(tfinal,'double') && isscalar(tfinal) ) %&& (tfinal>0) )
    error('tfinal must be a positive scalar number');
end

% Simulation Options
if ~isfield(opts,'Solver') || isempty(opts.Solver)
    % Use Local Events
    Solver = 'ode45';
else
    Solver = opts.Solver;
end
if ~isfield(opts,'Events') || isempty(opts.Events)
    % Use Local Events (This will wipe out Solver field)
    opts=odeset(opts,'Events',@LOCALevents);
end
opts.Solver = Solver;

if isempty(event_params)
    event_params.nbig = 1e6;
    event_params.nsmall = 1e-6;
    event_params.xbig = inf;
    event_params.xsmall = 0;
    event_params.funchandle = [];
end

% Pre-computation to avoid poly subsref during simulation.  Also
% need to re-order degmat based on order of variables in x vector
xv = char(x);
fdeg = f.degmat;
fcoef = f.coefficient;
fv = f.varname;
[tf,loc]=ismember(fv,xv);
fidx = loc;

% Run Simulations
xconv = zeros(N0,1);
xtraj = cell(N0,2);
for i1=1:N0
    x0 = x0mat(:,i1);
    event_params.nx0 = norm(x0);
    switch opts.Solver
        case 'ode15s'
            [t,x,te,xe,ie] = ...
                ode15s(@LOCALxdot,[0 tfinal],x0,opts,fdeg,fcoef,fidx,event_params);
        otherwise
            [t,x,te,xe,ie] = ...
                ode45(@LOCALxdot,[0 tfinal],x0,opts,fdeg,fcoef,fidx,event_params);
            
    end
    xtraj{i1,1} = t;
    xtraj{i1,2} = x;
    if ~isempty(ie) && (ie(1)==2 || ie(1)==4)
        xconv(i1)=1;
    end
end
xconv = logical(xconv);


%------------------------------------------------------------------------
% Local function which computes state derivative:
%       xdot = f(x)
% where f is a polynomial function of x.
%------------------------------------------------------------------------
function  xdot = LOCALxdot(t,x,fdeg,fcoef,fidx,event_params)

% Evaluate f using peval mex function
xdot = peval(x(fidx),fcoef,fdeg);
%if abs( t-floor(10*t)/10 ) < 1e-4, t, end

%------------------------------------------------------------------------
% Local event function
%    value(i) is the value of the ith event function.
%    isterminal(i) = 1 if the integration is to terminate at a zero of
%           this event function, otherwise, 0.
%    direction(i) = 0 if all zeros are to be located (the default),
%           +1 if only zeros where the event function is increasing,
%           and -1 if only zeros where the event function is decreasing.
%------------------------------------------------------------------------
function [value,isterminal,direction] = LOCALevents(t,x,fdeg,fcoef,fidx,event_params)

% Get event parameter data
nbig = event_params.nbig;
nsmall = event_params.nsmall;
xbig = event_params.xbig;
xsmall = event_params.xsmall;
funchandle = event_params.funchandle;
nx0 = event_params.nx0;

isterminal = ones(4,1);   % all events terminate the simulation
direction = zeros(4,1);   % terminate on all zeros of the events

% Simulation terminates if any value is equal to zero.
value = ones(4,1);
if isfinite(nbig)
    value(1) = (norm(x) <= nbig*nx0);
end
if nsmall>0
    value(2) = (norm(x) >= nsmall*nx0);
end
if isfinite(xbig)
    value(3) = all( abs(x) <= xbig );
end
if xsmall>0
    value(4) = any( abs(x) >= xsmall);
end

if ~isempty(funchandle)
    [value2,isterminal2,direction2]=funchandle(t,x,event_params);
    value = [value;value2];
    isterminal = [isterminal; isterminal2];
    direction = [direction; direction2];
end