function [sos,V] = sospolyvar(sos,ZSym,wscoeff)
% SOSPOLYVAR --- Declare a new polynomial variable in 
%       an SOS program 
%
% [SOSP,VAR] = sospolyvar(SOSP,Z)
%
% SOSP is the sum of squares program.
% VAR is the new polynomial variable.
% Z is the vector of monomials contained in VAR. Decision 
% variables corresponding to those monomials will be assigned 
% automatically by SOSPOLYVAR.
%
% Both VAR and Z are symbolic or polynomial objects.
%
% [SOSP,VAR] = sospolyvar(SOSP,Z,'wscoeff') will create
% the decision variables corresponding to VAR (i.e., coeff_xxx)
% also in MATLAB workspace.
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

if isfield(sos,'symvartable')

	if isnumeric(ZSym) & ZSym == 1
        ZSym = sym(ZSym);
	end;
	[sos,V] = sosvar(sos,'poly',ZSym);
	
	if nargin > 2 & wscoeff == 'wscoeff'
        var = sos.var.num;
        for i = sos.var.idx{var}:sos.var.idx{var+1}-1
            assignin('base',['coeff_',int2str(i-sos.var.idx{1}+1)],sos.symdecvartable(i));
        end;
	end;
    
else

    if isnumeric(ZSym) & ZSym==1
        pvar ZSym;
        ZSym = 0*ZSym+1;
    end;
    
    [dum,idx1,idx2] = intersect(ZSym.varname,sos.vartable);
    Z = sparse(size(ZSym.degmat,1),length(sos.vartable));
    Z(:,idx2) = sparse(ZSym.degmat(:,idx1));
    lenZ = size(Z,1);

    % Error handling needs to be added here, e.g. if Z is not a
    % valid monomial vector.
    
    % Add new variable
    sos.var.num = sos.var.num+1;
    var = sos.var.num;
    sos.var.type{var} = 'poly';
    sos.var.Z{var} = makesparse(Z);
    sos.var.ZZ{var} = makesparse(Z);
    sos.var.T{var} = speye(size(Z,1));
    sos.var.idx{var+1} = sos.var.idx{var}+size(Z,1);

    % Modify existing equations
    for i = 1:sos.expr.num
        sos.expr.At{i} = [sos.expr.At{i}; ...
            sparse(size(sos.var.T{var},1),size(sos.expr.At{i},2))];
    end;

    % Modify existing objective
    sos.objective = [sos.objective; sparse(sos.var.idx{var+1}-sos.var.idx{var},1)];        % 01/07/02

    % Modify decision variable table
    oldlen = length(sos.decvartable);
    sos.decvartable = [sos.decvartable; cell(lenZ,1)];
    for i = 1:lenZ
        sos.decvartable(oldlen+i) = {['coeff_',int2str(sos.var.idx{var}-sos.var.idx{1}+i)]};
    end;    
    
    pvar V;
    V = set(V,'varname',[sos.decvartable(oldlen+1:end); sos.vartable],...
        'degmat',[speye(lenZ) Z],...
        'coefficient',sparse(ones(lenZ,1)));
    
    if nargin > 2 & wscoeff == 'wscoeff'
        var = sos.var.num;
        for i = sos.var.idx{var}:sos.var.idx{var+1}-1
            pvar(['coeff_',int2str(i-sos.var.idx{1}+1)]);
            assignin('base',['coeff_',int2str(i-sos.var.idx{1}+1)],eval(['coeff_',int2str(i-sos.var.idx{1}+1)]));
        end;
	end;
    
    
end;
