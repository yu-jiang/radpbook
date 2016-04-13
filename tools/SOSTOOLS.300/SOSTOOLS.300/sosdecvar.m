function sos = sosdecvar(sos,decvartable)
% SOSDECVAR --- Declare new decision variables in an SOS 
%       program.
%
% SOSP = sosdecvar(SOSP,DECVARS)
%
% SOSP is the sum of squares program.
% DECVARS is a vector whose entries are the new decision 
% variables. DECVARS is either a symbolic or a polynomial object.
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


% 20/03/02 - SP

% Update tables and indices
if isfield(sos,'symvartable')

	decvartabletemp = sym2chartable(decvartable);  
	if length(sos.decvartable) ~= 2
        sos.decvartable = [decvartabletemp(1:end-1),',',sos.decvartable(2:end)];
        
    else
        sos.decvartable = [decvartabletemp(1:end)];  % Special case: previously empty table
	end;
	if size(decvartable,2) > 1
        decvartable = decvartable.';
	end;
	sos.symdecvartable = [decvartable; sos.symdecvartable];

else
    sos.decvartable = [decvartable.varname; sos.decvartable];
end;
    
    
offset = length(decvartable);
for i = 1:sos.var.num+1
    sos.var.idx{i} = sos.var.idx{i} + offset;
end;

% Update existing equations
for i = 1:sos.expr.num
    sos.expr.At{i} = [sparse(offset,size(sos.expr.At{i},2));sos.expr.At{i}]; 
end;

% Update existing objective
sos.objective = [sparse(offset,1); sos.objective];


