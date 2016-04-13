function Z3 = inconvhull(Z1,Z2)

% INCONVHULL --- Returns in Z3 the elements of Z1 that are
%                inside the convex hull of Z2
%
% This function returns the elements of Z1 that lie
% inside the convex hull of Z2. Z1 and Z2 must be matrices,
% each row of which is a point in n-space.
% Z1 is the set of monomials

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

% PP 20 Jan 03 % CDD
% AP 1 Feb 03 % convhulln + some other cases

% First, find the affine subspace where everything lives
% (for instance, in the homogeneous case)
nmons = size(Z2,1);

% Translate so it goes through the origin 
mr = mean(Z2);
Rzero = Z2 - repmat(mr,nmons,1) ;

% The columns of N generate the subspace 
N = null(Rzero);

% Z2*N should be constant 
cval = mean(Z2*N);

% Get only the monomials in the subspace
tol = .01 ;
ix = find(sum(abs(Z1*N-repmat(cval,size(Z1,1),1)),2)< tol) ;
nZ1 = Z1(ix,:);

% Now, the inequalities:
% Find an orthonormal basis for the subspace
% (I really should do both things at the same time...)

% Project to the lower dimensional space, so convhull works nicely

Q = orth(Rzero');

% This calls CDD, or whatever

if size(Z2*Q,2)>1
   if ~isempty(cddpath)
      [A,B] = vrep2hrep(Z2*Q);
   else 
      [A,B] = useconvhulln(Z2*Q);
   end
   % Find the ones that satisfy the inequalities, and keep them.
   ix = find(min(repmat(B,1,size(nZ1,1))-A*Q'*nZ1')>-tol);
   Z3 = nZ1(ix,:);
elseif size(Z2*Q,2)==1 % 1 Feb 03
   A = [1;-1];
   B = [max(Z2*Q);-min(Z2*Q)];
   ix = find(min(repmat(B,1,size(nZ1,1))-A*Q'*nZ1')>-tol);
   Z3 = nZ1(ix,:);
else
   Z3 = nZ1;
end

Z3 = unique(Z3,'rows'); % 1 Feb 03
