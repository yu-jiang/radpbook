function [Xsimdata] =  pplanesim(f,x,figno,x0,psimopts)
%function  [Xsimdata] = pplanesim(f,x,figno,x0,psimopts)
%
% DESCRIPTION
%   Draws the phase plane for a non-autonomous polynomial system:
%               dx/dt = f(x),   x(t) = x0
%
% INPUTS:
%   f: Vector field  (2-by-1 polynomial or function handle)
%   x: State  (2-by-1 vector of pvars
%   figno: Figure number for plotting
%   x0: 2-by-Npts array of initial conditions. Alternatively, the initial
%      conditions options can be specified as a structure with fields:
%      - range: 2-by-2 matrix with the i^th row specifying the min and
%        max value of the i^th state. default is [-1 1; -1 1]
%      - Npts: Number of initial conditions.  The actual number of points
%        depends on the sampling type (see sample below). (default is 100) 
%      - conv: True to plot only convergent trajectories (Default is false)
%      - div: True to plot only divergent trajectories (Default is false)
%      - sample: Sampling technique to be specified.  Choices are:
%         - 'grid': Generates ceil(sqrt(Npts)) points linearly spaced 
%            along each direction.
%         - 'bndry': Samples ceil(Npts/4) points along each of the
%            boundary specified by range.
%    psimopts:  Options structure passed to ODE solver. 
%
% OUTPUTS:
%    if no argument is invoked then only plot will be generated. However,
%    if one argument is invoked, then it will also return the simulation data.
%    For more information on the output refer to psim. The two outputs xtraj and
%    xconv will be bundled in the output argument as a cell array object.
%
% SYNTAX
%  pplanesim(f,x,figno,x0,psimopts)
%    Generate phase plane plot 
%  Xsimdata = pplanesim(f,x,figno,x0,psimopts)
%    Output all simulation data.


% Abhijit  12/01/09   Initial Coding


% XXX PJS 11/23/10
% 1) Why not have the same two output argument structure as psim?
% 2) Original help and code are written for f and x as generic Nx-by-1
%    These have to be 2-by-1 right?
% 3) Should figno be an option? In particular, I'm imagining the
%    structure of this function should be equivalent to psim:
%       [xtraj,xconv]=pplanesim(f,x,x0mat,tfinal,event_params,opts)
%    Except that opts now includes range, sample.   Is this the
%    right approach?  If so, then what is the point of also having
%    psim?  Let's discuss.  
% 4) The xlabel and ylabel should grab the varnames from x, i.e.
%    use labels x(1).varname and x(2).varname.
% 5) Spell out 'boundary' for the sample option.


% Error Checking

%---- only f can be provided only if f is a function handle
%---- f and x has to be provided if f is pvar.

if nargin == 1 && isa(f, 'function_handle')
    x = []; figno = 1;
elseif nargin == 1 && strcmp(class(f),'polynomial')
    error('x variables need to be provided if f is pvar object')
elseif nargin == 2 && strcmp(class(f),'polynomial')
    if ~isempty(setdiff(f.varname,char(x)))
        error('Variables do not match in f and x')
    else
        figno = 1;
    end
end
% --- set default for figno
if nargin == 2
    figno = 1;
end

% --------- Default x0
if nargin < 4
    x0 = [];
end

% --------- The user might not provide any psimopts.
if nargin < 5
    psimopts = [];
end

% -------- Handle x0

% --- Handle if x0 is data provided by user
if ~isstruct(x0) && ~isempty(x0)
    
    if (size(x0) ~= 2)
        error('incorrect dimension of data')
        % user provided npts -by-2 data
    elseif any(size(x0) == [ [] 2])== 1
        x0 = [x0(:,1)';x0(:,2)'];
    end
    
    % ---- check if its a structure or not and fill out information
elseif isstruct(x0) || isempty(x0)
    
    % --- Range
    if ~isfield(x0,'range') || isempty(x0.range)
        range = [-1 1; -1 1];
    else
        range = x0.range;
        % --------- range is not well defined
        if any((range(:,1) > range(:,2))== 1)
            error('Data range is not well defined')
        end
    end
    
    % --- Total number of points
    if ~isfield(x0,'Npts') || isempty(x0.Npts)
        Npts = 100;
        x0.Npts = Npts;
    else
        Npts = x0.Npts;
    end
    
    if ~isfield(x0,'conv') || isempty(x0.conv)
        plotconv = 0;
        x0.conv = plotconv;
    else
        plotconv = x0.conv;
    end
    
    if ~isfield(x0,'div') || isempty(x0.div)
        plotdiv = 0;
        x0.div = plotdiv;
    else
        plotdiv = x0.div;
    end
    
    if ~isfield(x0,'sample') || isempty(x0.sample)
        smpl = 'grid';
        x0.sample = 'grid';
    else
        smpl = x0.sample;
    end
    
end

% -------- Generate the IC if its not provided

xmin = range(1,1); xmax = range(1,2);
ymin = range(2,1); ymax = range(2,2);



if isstruct(x0)
    
    % ----- Sample Technique
    
    if strcmp(x0.sample,'grid')
        
        %---- MeshGrid data
        Ndpts = ceil(sqrt(Npts));
        x1 = linspace(xmin,xmax,Ndpts);
        y1 = linspace(ymin,ymax,Ndpts);
        [X,Y] = meshgrid(x1,y1);
        x0 = [X(:) Y(:)]';
        
    elseif strcmp(x0.sample,'bndry')
        
        Ndpts = ceil(Npts/4);
        x1 = linspace(xmin,xmax,Ndpts);
        y1 = linspace(ymin,ymax,Ndpts);
        X1Y1 = [min(x1)*ones(length(y1),1) y1'];
        X2Y2 = [max(x1)*ones(length(y1),1) y1'];
        X3Y3 = [x1' max(y1)*ones(length(x1),1)];
        X4Y4 = [x1' min(y1)*ones(length(x1),1)];
        x0 = [X1Y1; X2Y2; X3Y3; X4Y4]';
        
    end
    
end

% ------ Initialize storage data
if nargout ~= 0
    Xsimdata = cell(length(x0),2);
end

figure(figno)
Tfinal = 150; % Need to be incorporated either in psim or here

for j1 = 1:length(x0)
    
    % this calls ode45 for polynomial object
    if strcmp(class(f),'polynomial')
        [xtraj]=psim(f,x,x0(:,j1),Tfinal,psimopts);
        % Check if divergent, better method needed
        isdiv = norm(xtraj{2}(end,:)) >  5*norm(xtraj{2}(1,:));
        isconv = ~isdiv;
        
        % Plot results
        if isconv == 1 && plotdiv == 0
            plot(xtraj{2}(:,1), xtraj{2}(:,2),'-b'); hold on;
        elseif isdiv == 1 && plotconv == 0
            plot(xtraj{2}(:,1), xtraj{2}(:,2),'--r','LineWidth',2); hold on;
        end
        
    elseif isa(f, 'function_handle')
        [t,xtraj] = ode45(f,[0 Tfinal],x0(:,j1),psimopts);
        % Check if divergent, better method needed
        isdiv = norm(xtraj(end,:)) >  10*norm(xtraj(1,:));
        isconv = ~isdiv;
        
        % Plot results
        if isconv == 1 && plotdiv == 0
            plot(xtraj(:,1), xtraj(:,2),'-b'); hold on;
        elseif isdiv == 1 && plotconv == 0
            plot(xtraj(:,1), xtraj(:,2),'--r','LineWidth',2); hold on;
        end
    end
    
    axis([min(x0(1,:)) max(x0(1,:)) min(x0(2,:)) max(x0(2,:))]);
    drawnow
    
    if nargout ~= 0
        Xsimdata{j1,1} = xtraj;
        Xsimdata{j1,2} = isconv;
    end
    
    
end

% Create xlabel
xlabel('x','FontSize',14);
% Create ylabel
ylabel('y','FontSize',14);

