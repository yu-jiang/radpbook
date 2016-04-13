function [A,ZZ] = getconstraint(Z)
% GETCONSTRAINT --- Find constraint for sum of squares decomposition.
%
% [A,ZZ] = getconstraint(Z)
%
% Z is a monomial vector description.
% This function computes the constraint matrix A and the polynomial
% vector ZZ, such that if q satisfies
%
%    A*q = F'
%
% Then
%
%    Z'*Q*Z = F*ZZ      (where Q is the matrix form of q)
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

% 12/01/01 - SP

% Update: 12/12/01 -- Use sparse form. SP.
% Update: 24/12/01 -- Extra ZZ output. SP.

% First write Z'*Q*Z as
%
% Z'*Q*Z = (e1'*Q*R1 + e2'*Q*R2 + ... + en'*Q*Rn) ZZ

 
sizeZ = size(Z,1);
ZZ = Z + sprepmat(Z(1,:),sizeZ,1); 
M(1).R = speye(sizeZ);     % M contains the permutation matrices: R1, ..., Rn.
for i = 2:sizeZ
    Ztemp = Z + sprepmat(Z(i,:),sizeZ,1);
    [R1,R2,ZZ] = findcommonZ(ZZ,Ztemp);
    for j = 1:i-1
        M(j).R = M(j).R*R1;
    end;
    M(i).R = R2;
end;
 

% Construct the constraint equations
Q = sparse([],[],[],sizeZ,sizeZ,1);
A = sparse(size(ZZ,1),sizeZ^2);         % 12/12/01
for i = 1:(sizeZ)^2
    Q(i) = 1; 
    %for j = 1:sizeZ
    [j,k] = find(Q);                    % 12/12/01 
    A(:,i) = A(:,i) + M(j).R' * Q(j,:)';
    %end;
    Q(i) = 0;
end;
 