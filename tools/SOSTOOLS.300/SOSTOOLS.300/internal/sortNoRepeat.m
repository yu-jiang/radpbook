function [Znew] = sortNoRepeat(Z,Zin)
%This function excludes rows of Zin which are already in matrix Z
%Z and Zin have the same dimension
% Z,Zin do not contain repeated rows 

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

if isempty(Zin)
    Znew = Z; return
end
if isempty(Z)
    Znew = flipud(sortrows([Zin])); return
end


[nmonZ,nv] = size(Z);

[nmonZin,ncolzin] = size(Zin);
Znew = [];
for k = 1:nmonZin

    check = sum(ismember(Z,Zin(k,:),'rows'));
    if check==0
        Znew = [Znew;Zin(k,:) ];
    end
end

Znew = flipud(sortrows([Z; Znew]));