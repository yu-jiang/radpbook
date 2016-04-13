function sos = sossetobj(sos,symexpr)
% SOSSETOBJ --- Set the objective function of an SOS program. 
%
% SOSP = sossetobj(SOSP,EXPR)
%
% Given a sum of squares program SOSP and an expression EXPR, SOSSETOBJ
% makes EXPR the objective function of the sum of squares program SOSP, 
% i.e., EXPR will be minimized.
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


% 01/07/02 - SP
% 03/01/02 - SP -- New syntax

if isfield(sos,'symvartable')

	% Creating matching table
	objvartable = [',',findsym(symexpr),','];
	objvartable = objvartable(find(objvartable~=' '));        % 03/25/02
	decvartable = [',',sos.decvartable(2:end-1),','];
	idxcomma1 = find(objvartable==',');
	idxcomma2 = find(decvartable==',');
	
	for i = 1:length(idxcomma1)-1
        j = findstr(decvartable,objvartable(idxcomma1(i):idxcomma1(i+1)));
        if length(j)==0
            error('There is unknown decision variable');
        end;
        j = find(idxcomma2==j);
        idx(i) = j;        % idx contains indices of decision variables
        tempexpr = diff(symexpr,sos.symdecvartable(idx(i)));
        if ~isempty(findsym(tempexpr))
            error('Not a linear objective function.');
        else
            sos.objective(idx(i)) = double(tempexpr);
        end;
        symexpr = simplify(symexpr - sos.objective(idx(i))*sos.symdecvartable(idx(i)));    % Optimize here
	end;
    
else
    [dummy,idxdecvar1,idxdecvar2] = intersect(symexpr.varname,sos.decvartable);
    for i = 1:length(idxdecvar2)
        sos.objective(idxdecvar2(i)) = symexpr.coefficient(idxdecvar1(i));
    end;
    
    
end;