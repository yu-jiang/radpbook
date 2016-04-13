function [GAM,vars,xopt] = findbound(p,ineq,eq,DEG,options)
% FINDBOUND --- Find a global/constrained lower bound for a polynomial. 
%
% [GAM,VARS,XOPT] = findbound(P,OPTIONS)
%
% P is a polynomial of even degree, whose global lower bound is to be
% computed. The function computes the largest GAM such that (P - GAM) is a
% sum of squares.
%
% OPTIONS is an optional argument which specifies the solver and the
% solver-specific parameters on its fields as
%   options.solver
%   options.params
% the default value for options.solver is 'SeDuMi'. The solver parameters
% are fields of options.params as, for instance, options.params.tol = 1e-9.
%
%
% In addition, a vector VARS containing all the variables in P, and the
% optimal argument XOPT (in the same ordering as VARS) will  be computed,
% such that by substituting VARS = XOPT into P  we obtain P = GAM. If no
% such optimal argument is found, the function  will return an empty XOPT.
%
% [GAM,VARS,XOPT] = findbound(P,INEQ,EQ,DEG,OPTIONS) will compute a lower bound 
% for the constrained optimization problem: minimize P subject to INEQ >= 0
% and EQ = 0. Here INEQ and EQ are cells of polynomials which define the
% inequality and equality constraints. The function computes a lower bound
% based on Schmudgen's positivstellensatz, i.e., it computes the largest
% GAM such that
% 
%   (P-GAM) - LAM'*EQ - SIGMA1'*INEQ - ....
% 
% is a sum of squares, where LAM is a vector of polynomials, and SIGMA is a
% vector of sums of squares. The degree of the expression will be
% determined by the input argument DEG. Higher value of DEG will return a
% better lower bound, although the computational cost will also increase.
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

switch nargin
    case 1 
        options.solver='sedumi';
        ineq = [];
        eq = [];
    case 2 
        options = ineq;
        ineq = [];
        eq = [];
    case 3
        options.solver='sedumi';
        ineq = ineq(:);
        DEG = eq;
        eq = [];
    case 4
        if ~isstruct(DEG)
            ineq = ineq(:);
            eq = eq(:);
            options.solver='sedumi';
        else
            options = DEG;
            ineq = ineq(:);
            DEG=eq;
            eq = [];
        end
    case 5 
        ineq = ineq(:);
        eq = eq(:);
end


vect = [p; ineq; eq];

% Find the independent variables, check the degree
if isa(vect,'sym')
   varschar = findsym(vect);
   vars = sym(['[',varschar,']']);
   nvars = size(vars,2) ; 
   if nargin > 2
       degree = 2*floor(DEG/2);
       deg = zeros(length(vect),1);
       for i = 1:length(vect)
           deg(i) = double(feval(symengine,'degree',vect(i),converttochar(vars)));
           if deg(i) > degree
               error('One of the expressions has degree greater than DEG.');
           end;
       end;   
   else
       % Can change, to call maple only once.
       for var = 1:nvars;
           if rem(double(feval(symengine,'degree',p)),2) ;
               disp(['Degree in ' char(vars(var)) ' should be even. Otherwise the polynomial is unbounded.']);
               GAM = -Inf;
               xopt = [];
               return;
           end;
       end;
   end;
   syms gam;
   
else
   varname = vect.var;
   vars = [];
   for i = 1:size(varname,1)
       pvar(varname{i});
       vars = [vars eval(varname{i})];
   end;
   
   if nargin > 2
       degree = 2*floor(DEG/2);
       degmat = sum(vect.degmat,2);
       deg = zeros(length(vect),1);
       for i = 1:length(vect)
           idx = find(vect.coefficient(:,i));
           deg(i) = max(degmat(idx));
           if deg(i) > degree
               error('One of the expressions has degree greater than DEG.');
           end;
       end;   
   else
       deg = mod(max(p.degmat,[],1),2);
       if sum(deg)~=0
           i = find(deg~=0);
           disp(['Degree in ' varname{i(1)} ' should be even. Otherwise the polynomial is unbounded.']);
           GAM = -Inf;
           xopt = [];
           return;
       end;
   end;
   pvar gam;
end;

% Construct other valid inequalities
if length(ineq)>1
    for i = 1:2^length(ineq)-1
        Ttemp = dec2bin(i,length(ineq));
        T(i,:) = str2num(Ttemp(:))';
    end;
    T = T(find(T*deg(2:1+length(ineq)) <= degree),:);
    
    deg = [deg(1); T*deg(2:1+length(ineq)); deg(2+length(ineq):end)];
    for i = 1:size(T,1)
        ineqtempvect = (ineq.').^T(i,:);
        ineqtemp(i) = ineqtempvect(1);
        for j = 2:length(ineqtempvect)
            ineqtemp(i) = ineqtemp(i)*ineqtempvect(j);
        end;
    end;
    ineq = ineqtemp;
end;

prog = sosprogram(vars,gam);
expr = p-gam;
for i = 1:length(ineq)
    [prog,sos] = sossosvar(prog,monomials(vars,0:floor((degree-deg(i+1))/2)));
    expr = expr - sos*ineq(i);
end;
for i = 1:length(eq)
    [prog,pol] = sospolyvar(prog,monomials(vars,0:degree-deg(i+1+length(ineq))));
    expr = expr - pol*eq(i);
end;
prog = sosineq(prog,expr);
prog = sossetobj(prog,-gam);
[prog,info] = sossolve(prog,options);


xopt = [];

if (info.dinf>1e-2)|(info.pinf>1e-2)
    disp('No lower bound could be computed. Unbounded below or infeasible?');
    GAM = '-inf';
    return;
else
    GAM = double(sosgetsol(prog,gam,16));
    % Returning the solution
    % The code here will only work in the rank one case.
    % Otherwise, we'll solve the f_i(x)=0 (efficiently)
    %

    % Find where the variables are (perhaps they're in the wrong order? Ask Stephen.
    ix = find(sum(prog.extravar.Z{1},2)==1) ;
    xopt = prog.solinfo.extravar.dual{1}(ix,1) ;
    vars = vars.';
    
    % If the upper and lower bounds are close (absolute or relative), return them
    ach = double(subs(p,num2cell(vars).',num2cell(xopt).'));
    
    if min(abs(ach/GAM-1),abs(ach-GAM)) > 1e-4 ; 
        xopt = [];
        return;
    end
    
    % Check inequality and equality constraints
    if length(ineq)>1
        ineq = double(subs(ineq,num2cell(vars).',num2cell(xopt).'));
    else
        ineq = 0;
    end;
    if length(eq)>1
        eq = double(subs(eq,num2cell(vars).',num2cell(xopt).'));
    else
        eq = 0;
    end;
    if min(ineq)<-1e-6 | max(abs(eq))>1e-6
        xopt = [];
        return;
    end;
    

    
end;



