function Z = sparsemultipart(Z1,Z2,info)

% Find the elements in Z1 that are in the convex hull of Z2, where Z2 is bipartite

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

% AP Feb 01 0

[Y,I] = find(sum(Z2,1)==0);
sizeinf = length(info);
  
if sizeinf == 1
   error('Error in sparsemultipart option - at least two sets are required');
end 

for i = 1:sizeinf
  Z3{i} = inconvhull(Z1(:,info{i}),Z2(:,info{i}));
end

for i = 1:sizeinf-1%GVcomment for the variables introducing homogeneous monomials 
  Z3{i+1} = [repmat(Z3{i},size(Z3{i+1},1),1),kron(Z3{i+1},ones(size(Z3{i},1),1))];
end
Z3 = Z3{sizeinf};
Z = zeros(size(Z3,1),size(Z2,2));
lgth = 0;
for i = 1:sizeinf
   Z(:,info{i})= Z3(:,(1+lgth):(lgth+length(info{i})));
   lgth = length(info{i})+lgth;
end
Z(:,I) = zeros(size(Z,1),size(I,2));