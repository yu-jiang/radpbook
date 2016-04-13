function [sos,V] = sosvar(sos,Type,ZSym)
% SOSVAR --- Declare a new SOSP variable in an SOS program 
%
% [SOSP,VAR] = sosvar(SOSP,TYPE,Z)
%
% SOSP is the sum of squares program.
% VAR is the new SOSP variable. This variable is described by TYPE and a 
% vector of monomial Z in the following manner:
%   TYPE = 'sos'  : the new SOSP variable is a sum of squares        --- Z' * Q * Z
%   TYPE = 'poly' : the new SOSP variable is an ordinary polynomial  --- F * Z
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


% 12/24/01 - SP
% 01/07/02 - SP
% 02/21/02 - SP -- Symbolic polynomial
% 03/01/02 - SP -- Symbolic V

if (~strcmp(Type,'sos') & ~strcmp(Type,'poly'))
    error('Unknown type.');
end;

Z = getdegrees(sym2chartable(ZSym),sos.vartable);

% Error handling needs to be added here, e.g. if Z is not a
% valid monomial vector.
    
% Add new variable
sos.var.num = sos.var.num+1;
var = sos.var.num;
sos.var.type{var} = Type;
sos.var.Z{var} = makesparse(Z);
switch Type
case 'sos'
    [T,ZZ] = getconstraint(Z);    
    sos.var.ZZ{var} = ZZ;
    sos.var.T{var} = T';
    sos.var.idx{var+1} = sos.var.idx{var}+size(Z,1)^2;
case 'poly'
    sos.var.ZZ{var} = makesparse(Z);
    sos.var.T{var} = speye(size(Z,1));
    sos.var.idx{var+1} = sos.var.idx{var}+size(Z,1);
end;

% Modify existing equations
for i = 1:sos.expr.num
    sos.expr.At{i} = [sos.expr.At{i}; ...
            sparse(size(sos.var.T{var},1),size(sos.expr.At{i},2))];
end;

% Modify existing objective
sos.objective = [sos.objective; sparse(sos.var.idx{var+1}-sos.var.idx{var},1)];        % 01/07/02

% Modify decision variable table
decvartabletemp = [];
for i = sos.var.idx{var}:sos.var.idx{var+1}-1
    decvartabletemp = [decvartabletemp,'coeff_',int2str(i-sos.var.idx{1}+1),','];
end;
decvartabletemp = decvartabletemp(1:end-1);
if length(sos.decvartable)==2
    sos.decvartable = sos.decvartable(1);
end;
sos.decvartable = [strrep(sos.decvartable,']',','),decvartabletemp,']'];
sos.symdecvartable = [sos.symdecvartable; ...
        sym(['[',strrep(decvartabletemp,',',';'),']'])];

% Return variable expression -- Possibly faster using Maple
switch Type
case 'sos'
    if size(sos.var.ZZ{var},1) <= 10
        % This is slow, but give clear visualization.
        ZZSym = mysympower(sos.symvartable,sos.var.ZZ{var});
        coeff = myAdecvarprod(sos.var.T{var},sos.symdecvartable(sos.var.idx{var}:sos.var.idx{var+1}-1));
        V = (coeff.')*ZZSym;
    else
        % This is fast
        dummy = sos.symdecvartable(sos.var.idx{var}:sos.var.idx{var+1}-1);
        V = ZSym.' * reshape(dummy,sqrt(length(dummy)),sqrt(length(dummy))) * ZSym;
    end;
case 'poly'
    V = sos.symdecvartable(sos.var.idx{var}:sos.var.idx{var+1}-1).' * ZSym;
end;

%on = maple('convert',1,'float');   % 08/07/02
%V = on*V;                          % 08/07/02

% ===============================================
function sympoly = myAdecvarprod(A,decvar);

charpoly = 'matrix([';
for i = 1:size(A,2)
    idx = find(A(:,i));
    if ~isempty(idx)
        charpolytemp = [];
        for j = idx'
            charpolytemp = [charpolytemp,char(decvar(j)),'+'];
        end;
        charpoly = [charpoly,'[',charpolytemp(1:end-1),'], '];
    else
        charpoly = [charpoly,'[0], '];
    end;
end;
charpoly = [charpoly(1:end-2),'])'];
sympoly = sym(charpoly);

