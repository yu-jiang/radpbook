function V = findlyap(f,vars,degree,options)
% FINDLYAP --- Find a Lyapunov function V(x) for a simple polynomial 
%        dynamical system.
%
% V = findlyap(F,VARS,DEGREE,OPTIONS)
%
% F is the (polynomial) vector field  of the dynamical system. VARS is the
% variables of the system. F and VARS must be ordered such that
%
%    dVARS/dt = F
%
% DEGREE is the desired degree of the Lyapunov function (has to be  an even
% integer).
%
% OPTIONS is an optional argument which specifies the solver and the
% solver-specific parameters on its fields as
%   options.solver
%   options.params
% the default value for options.solver is 'SeDuMi'. The solver parameters
% are fields of options.params as, for instance, options.params.tol = 1e-9.
%
%
% A polynomial Lyapunov function V will be computed such that
% 
%    V is positive definite, radially unbounded. dV/dt is negative
%    semidefinite along the system trajectory.
% 
% If no such Lyapunov function exists, the function will return  empty V.
%

% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 3.00.
%
% Copyright (C)2002, 2004, 2013  A. Papachristodoulou (1), J. Anderson (1),
%                                G. Valmorbida (1), S. Prajna (2), 
%                                P. Seiler (3), P. A. Parrilo (4)
% (1) Department of Engineering Science, University of Oxford, Oxford, U.K.
% (2) Control and Dynamical Systems - California Institute of Technology,
%     Pasadena, CA 91125, USA.
% (3) Aerospace and Engineering Mechanics Department, University of
%     Minnesota, Minneapolis, MN 55455-0153, USA.
% (4) Laboratory for Information and Decision Systems, M.I.T.,
%     Massachusetts, MA 02139-4307
%
% Send bug reports and feedback to: sostools@cds.caltech.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%


% 12/27/01 - SP
% 03/27/02 - SP

if nargin<4
    %options.solver = 'sedumi';
    options.solver = 'sdpt3';
end

if degree<2
    error('Order has to be at least 2.');
end;
if rem(degree,2)
    error('Order is not even.');
end;
if length(vars)~=length(f)
    error('F and VARS have different length.');
end;

if size(vars,1)>1
    vars = vars.';
end;

prog = sosprogram(vars);
[prog,V] = sospolyvar(prog,monomials(vars,[2:degree]));

% Positive definiteness condition
expr1 = V;
for i = 1:length(vars)
    expr2 = -0.1;
    for j = 2:2:degree
        [prog,epsmat(i,j/2)] = sossosvar(prog,1);
        expr2 = expr2 + epsmat(i,j/2);
    end;
    prog = sosineq(prog,expr2);
    if isfield(prog,'symvartable')
       expr1 = expr1 - epsmat(i,:) * mysympower(vars(i),[2:2:degree]');
    else
       expr1 = expr1 - epsmat(i,:) * (vars(i).^[2:2:degree]');
    end;
end;
prog = sosineq(prog,expr1);

% Condition on derivative
expr1 = 0;
for i = 1:length(vars)
    expr1 = expr1 - diff(V,vars(i))*f(i);
end;
prog = sosineq(prog,expr1);

[prog,info] = sossolve(prog,options);
if (info.dinf==1)|(info.pinf==1)
    disp('No Lyapunov function has been found for this system.');
    V = [];
else
    V = sosgetsol(prog,V);
end;





