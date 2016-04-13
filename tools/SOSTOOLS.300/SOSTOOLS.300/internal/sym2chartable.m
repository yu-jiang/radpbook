function chartable = sym2chartable(symtable)
% SYM2CHARTABLE -- Convert symbolic vector to string table.
% 
% STR = sym2chartable(SYM)
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

% 03/01/02 - SP

chartable = char(symtable);

if length(symtable) == 1
    verstring = char(version);
    if verstring(1) =='7' || verstring(1) =='8'
        chartable = ['matrix(',chartable,')'];
    else
        chartable = ['array(',chartable,')'];
    end;
end;
chartable = strrep(chartable,'[','');
chartable = strrep(chartable,']','');
verstring = char(version);
if verstring(1) =='7' || verstring(1) =='8'
    chartable = strrep(chartable,'matrix(','[');
else
    chartable = strrep(chartable,'array(','[');
end;
chartable = strrep(chartable,')',']');
chartable = strrep(chartable,' ','');  %JA 03/2013: for compatibility with R2011b
