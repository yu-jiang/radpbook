function [sos,V] = sossosvar(sos,ZSym,wscoeff)
% SOSSOSVAR --- Declare a new sum of squares variable in 
%       an SOS program 
%
% [SOSP,VAR] = sossosvar(SOSP,Z)
%
% SOSP is the sum of squares program.
% VAR is the new sum of squares variable.
% ZSym is the column vector of monomials contained in the Gram 
% decomposition of VAR, i.e.,
%
%       VAR = ZSym' * COEFF * ZSym
%
% where COEFF is a coefficient matrix that is restricted to be 
% positive semidefinite. COEFF and the decision variables contained 
% in it will be constructed automatically by SOSSOSVAR.
%
% Both VAR and Z are either symbolic or polynomial objects.
%
% [SOSP,VAR] = sospolyvar(SOSP,ZSym,'wscoeff') will create
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

if isfield(sos,'symvartable')

	if isnumeric(ZSym) & ZSym == 1
        ZSym = sym(ZSym);
	end;
	[sos,V] = sosvar(sos,'sos',ZSym);
	
	if nargin > 2 & wscoeff == 'wscoeff'
        var = sos.var.num;
        for i = sos.var.idx{var}:sos.var.idx{var+1}-1
            assignin('base',['coeff_',int2str(i-sos.var.idx{1}+1)],sos.symdecvartable(i));
        end;
	end;
    
else
    
    
    if isnumeric(ZSym) & ZSym == 1
	pvar ZSym;
        ZSym = ZSym*0+1;
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
    sos.var.type{var} = 'sos';
    sos.var.Z{var} = makesparse(Z);
    [T,ZZ] = getconstraint(Z);    
    sos.var.ZZ{var} = ZZ;
    sos.var.T{var} = T';
    sos.var.idx{var+1} = sos.var.idx{var}+size(Z,1)^2;

    % Modify existing equations
    for i = 1:sos.expr.num
        sos.expr.At{i} = [sos.expr.At{i}; ...
            sparse(size(sos.var.T{var},1),size(sos.expr.At{i},2))];
    end;

    % Modify existing objective
    sos.objective = [sos.objective; sparse(sos.var.idx{var+1}-sos.var.idx{var},1)];        % 01/07/02

    % Modify decision variable table
    oldlen = length(sos.decvartable);
    sos.decvartable = [sos.decvartable; cell(lenZ^2,1)];
    for i = 1:lenZ^2
        sos.decvartable(oldlen+i) = {['coeff_',int2str(sos.var.idx{var}-sos.var.idx{1}+i)]};
    end;    
    
    pvar V;

    ZTemp = sparse(lenZ^2,size(Z,2));
    for i = 1:lenZ
        ZTemp((i-1)*lenZ+[1:lenZ],:) = Z + sprepmat(Z(i,:),lenZ,1);
    end;
    V = set(V,'varname',[sos.decvartable(oldlen+1:end); sos.vartable],...
        'degmat',[speye(lenZ^2) ZTemp],...
        'coefficient',sparse(ones(lenZ^2,1))); % MMP 6.19.2013 - fixed compatability bug with POLVAR 2013
    
    if nargin > 2 & wscoeff == 'wscoeff'
        var = sos.var.num;
        for i = sos.var.idx{var}:sos.var.idx{var+1}-1
            pvar(['coeff_',int2str(i-sos.var.idx{1}+1)]);
            assignin('base',['coeff_',int2str(i-sos.var.idx{1}+1)],eval(['coeff_',int2str(i-sos.var.idx{1}+1)]));
        end;
	end;
end;
