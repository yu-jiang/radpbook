function p = sosgetsol(sos,V,digit)
% SOSGETSOL --- Get the solution from a solved SOS program 
%
% SOL = sosgetsol(SOSP,VAR,DIGIT) 
%
% SOL is the solution from a solved sum of squares program SOSP,
% obtained through substituting all the decision variables
% in VAR by the numerical values which are the solutions to
% the corresponding semidefinite program. 
%
% The third argument DIGIT (optional) will determine the 
% accuracy of SOL in terms of the number of digits. Default 
% value is 5.
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

% 12/27/01 - SP
% 02/21/02 - SP -- Symbolic polynomial
% 03/01/02 - SP -- New syntax
% 03/15/02 - SP -- Fast

if nargin == 2
    digit = 5;   % Default
end;

if isfield(sos,'symvartable')
    
    p = mysymsubs(V,sos.symdecvartable,sos.solinfo.RRx(1:length(sos.symdecvartable)),digit);

else
    
    [dummy,idxdecvar1,idxdecvar2] = intersect(V.varname,sos.decvartable);
    idxvar = setdiff(1:length(V.varname),idxdecvar1);
    coeffs = V.degmat(:,idxdecvar1)*sos.solinfo.RRx(idxdecvar2);
    coeffs = V.coefficient.*coeffs;         % 01/31/02
    varname = V.varname(idxvar);
    if isempty(idxvar)
        degmat = 0;
    else
        degmat = V.degmat(:,idxvar);
    end;
    if isempty(varname) 
        varname = cellstr('dummyvar');
    end
    p = set(V,'varname',varname,'degmat',degmat,'coefficient',coeffs);
end;
