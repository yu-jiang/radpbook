function [Q,Z,decomp,Den] = findsos(P,flag,options)
% FINDMATRIXSOS --- Find a sum of squares decomposition of a given matrix polynomial.
%
% [Q,Z,decomp,Den] = findsos(P,flag,options)
%
% P is a symmetric polynomial matrix.
%
% FLAG is an optional argument which, if set to 'rational', returns a
% sum of squares decomposition with rational coefficients. 
%
% OPTIONS is an optional argument which specifies the solver and the
% solver-specific parameters on its fields as
%   options.solver
%   options.params
% the default value for options.solver is 'SeDuMi'. The solver parameters
% are fields of options.params as, for instance, options.params.tol = 1e-9.
%
% A positive semidefinite Q and a symbolic monomial vector Z will be
% computed such that
%
%    (Ir kron Z)' * Q * (Ir kron Z) = P(x)
%
% If P is not a sum of squares, the function will return empty Q and Z.
%
% If P is a polynomial with integer coefficients and is represented as a
% symbolic object, then [Q,Z] = findsos(P,'rational') will compute a rational
% matrix Q such that Z'*Q*Z = P.

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


if nargin == 2 
    if ~strcmp(flag,'rational')
        options=flag;
        flag = 'abcdefgh';
    else
        %options.solver='sedumi';
        options.solver='sdpt3';
        flag = 'abcdefgh';
    end
elseif nargin==1
    %options.solver='sedumi';
    options.solver='sdpt3';
end

if isa(P,'sym')
    
    dimp = size(P);
    if dimp(1)~=dimp(2) ;
        disp(['The polynomial matrix is not square, it cannot be a sum of squares']);
        Q = [];
        Z = [];
        return;
    end;
    
    
    P = expand(P);
    vars = findsym(P);
    vars = sym(['[',vars,']']);
    
    nvars = size(vars,2) ;
    
    
    for var = 1:nvars;
        for j = 1:dimp(1)
            if rem(double(feval(symengine,'degree',P(j,j))),2) ;
                disp(['Degree in ' char(vars(var)) ' is not even for the diagonal elements. The polynomial matrix cannot be a sum of squares']);
                Q = [];
                Z = [];
                return;
            end;
        end;
    end;
    
    if isequal(P,P.')~=1;
        disp(['The polynomial matrix is not symmetric, it can not be a sum of squares']);
        Q = [];
        Z = [];
        return;
    end;
    
    
    
    prog = sosprogram(vars);
    prog = sosineq(prog,P);
    [prog,info] = sossolve(prog,options);
    
    if strcmp(options.solver,'SDPNAL')
        disp('findsos function currently not supported for SDPNAL solver');
        Q = [];
        Z = [];
        if nargout >= 3
            decomp = [];
            Den = [];
        end
        return;
    elseif (info.dinf==1)||(info.pinf==1)
        disp('No sum of squares decomposition is found.');
        Q = [];
        Z = [];
        if nargout >= 3
            decomp = [];
            Den = [];
        end
        return;
    else
        Q = reshape(prog.solinfo.RRx,sqrt(length(prog.solinfo.RRx)),sqrt(length(prog.solinfo.RRx)));
        if isa(P,'sym');
            Z = mysympower(vars,prog.extravar.Z{1});
        else
            pvar Z
            coefftemp = speye(size(prog.extravar.Z{1},1));
            Z = polynomial(coefftemp,prog.extravar.Z{1},prog.vartable,[size(prog.extravar.Z{1},1),1]);
        end;
    end;
    
    if nargin > 1 & strcmp(flag,'rational')
        
        A = (full(prog.expr.At{1})') ;
        B = (full(prog.expr.b{1})) ;
        % Round x0
        N = 1 ; % This is the tolerance in the rounding
        xx = prog.solinfo.x;
        
        kmax = 10 ;
        pd = 0 ;
        while kmax ~= 0 ;
            kmax = kmax - 1 ;
            x0 = round(N*xx);
            [Qflat,NN] = proj3(x0,A,B,N);
            n = sqrt(length(Qflat));
            Qr = reshape(Qflat,n,n);
            % Qr should be PSD (should really check symbolically)
            if min(eig(Qr/NN))>-1e-14 ; kmax=0 ; pd = 1 ; end
            % Increase N, and try again
            N = 2*N ;
        end
        
        % Experimental, no good error checking yet, so we check that
        % expand(NN*P - Z.'*Qr*Z) is zero!
        
        if (expand(NN*P-Z.'*Qr*Z) ~= 0) | (pd == 0);
            Qr=[];Z=[];NN=[];
            disp('Could not compute a rational SOS!');
        end
        
        if nargout == 4
            Q = Qr;
            Den = NN;
        else
            if isa(P,'sym')
                Q = sym(Qr/NN);
            else
                Q = Qr/NN;
                disp('To obtain an exact rational solution, run the function with three output arguments.');
            end;
        end;
        
    end;
    
    
    L = sqrtm(double(Q));
    decomp   = L*(kron(eye(dimp(1)),Z));

else 

    
    % PJS -- Handle polynomial variables
    % Most of this code simply mimics the symbolic case and hence the
    % overlap can be reduced.

    dimp = size(P);
    if dimp(1)~=dimp(2) ;
        disp(['The polynomial matrix is not square, it cannot be a sum of squares']);
        Q = [];
        Z = [];
        return;
    end;

    nvars = P.nvar;
    vars = polynomial(zeros(nvars,1));
    for j=1:nvars
        vars(j) = pvar(P.varname{j});
    end
    
    for var = 1:nvars;
        for j = 1:dimp(1)
            Pjj = P(j,j);
            maxdeg = max(Pjj.degmat(:,var));
            if rem(maxdeg,2)
                disp(['Degree in ' char(vars(var)) ' is not even for the diagonal elements. The polynomial matrix cannot be a sum of squares']);
                Q = [];
                Z = [];
                return;
            end;
        end;
    end;

    tmp = isequal(P,P.');
    if ~all(tmp(:))
        disp(['The polynomial matrix is not symmetric, it can not be a sum of squares']);
        Q = [];
        Z = [];
        return;
    end;
        
    prog = sosprogram(vars);
    prog = sosineq(prog,P);
    [prog,info] = sossolve(prog);
    
    if 'solver' == 'SDPNAL'
        disp('findsos function currently not supported for SDPNAL solver');
        Q = [];
        Z = [];
        if nargout >= 3
            decomp = [];
            Den = [];
        end
        return;
    elseif (info.dinf==1)||(info.pinf==1)
        disp('No sum of squares decomposition is found.');
        Q = [];
        Z = [];
        if nargout >= 3
            decomp = [];
            Den = [];
        end
        return;
    else
        Q = reshape(prog.solinfo.RRx,sqrt(length(prog.solinfo.RRx)),sqrt(length(prog.solinfo.RRx)));
        if isa(P,'sym');
            Z = mysympower(vars,prog.extravar.Z{1});
        else
            pvar Z
            coefftemp = speye(size(prog.extravar.Z{1},1));
            Z = polynomial(coefftemp,prog.extravar.Z{1},prog.vartable,[size(prog.extravar.Z{1},1),1]);
        end;
    end;
    
    if nargin > 1 & strcmp(flag,'rational')
        
        A = (full(prog.expr.At{1})') ;
        B = (full(prog.expr.b{1})) ;
        % Round x0
        N = 1 ; % This is the tolerance in the rounding
        xx = prog.solinfo.x;
        
        kmax = 10 ;
        pd = 0 ;
        while kmax ~= 0 ;
            kmax = kmax - 1 ;
            x0 = round(N*xx);
            [Qflat,NN] = proj3(x0,A,B,N);
            n = sqrt(length(Qflat));
            Qr = reshape(Qflat,n,n);
            % Qr should be PSD (should really check symbolically)
            if min(eig(Qr/NN))>-1e-14 ; kmax=0 ; pd = 1 ; end
            % Increase N, and try again
            N = 2*N ;
        end
        
        % Experimental, no good error checking yet, so we check that
        % expand(NN*P - Z.'*Qr*Z) is zero!
        
        if (expand(NN*P-Z.'*Qr*Z) ~= 0) | (pd == 0);
            Qr=[];Z=[];NN=[];
            disp('Could not compute a rational SOS!');
        end
        
        if nargout == 4
            Q = Qr;
            Den = NN;
        else
            if isa(P,'sym')
                Q = sym(Qr/NN);
            else
                Q = Qr/NN;
                disp('To obtain an exact rational solution, run the function with three output arguments.');
            end;
        end;
        
    end;
    
    
    L = sqrtm(double(Q));
    decomp   = L*(kron(eye(dimp(1)),Z));
    
    
    % PJS -- Comment out original code in this section
%     varname = P.var;
%     vars = [];
%     for i = 1:size(varname,1)
%         pvar(varname{i});
%         vars = [vars eval(varname{i})];
%     end;
%     
%     deg = mod(max(P.degmat,[],1),2);
%     if sum(deg)~=0
%         i = find(deg~=0);
%         disp(['Degree in ' varname{i(1)} ' is not even. The polynomial cannot be a sum of squares']);
%         Q = [];
%         Z = [];
%         return;
%     end;
    
end;