function [R1,R2,Z] = findcommonZ(Z1,Z2)
% FINDCOMMONZ --- Find common Z and permutation matrices R1, R2
%
% [R1,R2,Z] = findcommonZ(Z1,Z2)
%
% Given two vectors of monomials Z1 and Z2, this 
% function will compute another vector of monomials Z
% containing all the monomials of Z1 and Z2, and
% permutation matrices R1, R2 such that
%
%  Z1 = R1*Z
%  Z2 = R2*Z
%
% Assumption: all the monomials in Z1, as well as
% the monomials in Z2, are DISTINCT --- but Z1 and Z2 may 
% have common monomials.
%
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

% 11/26/01 - SP

% Update: 12/12/01 -- Use sparse form
% Update: 12/26/01 -- Empty and identical monomials

if (size(Z1,1) + size(Z2,1)) <= 1           % 12/26/01 -- Empty monomials
    Z = [Z1; Z2];
    R1 = speye(size(Z1,1),size(Z,1));
    R2 = speye(size(Z2,1),size(Z,1));
    return;
end;

% Constructing Index matrix
sizeZ1 = size(Z1,1);
Ind1 = (1:sizeZ1)';
sizeZ2 = size(Z2,1);
Ind2 = (1:sizeZ2)';
Ind = [Ind1 zeros(size(Ind1))+sizeZ2+1; zeros(size(Ind2))+sizeZ1+1 Ind2];

% Constructing Z
ZZ = [Z1; Z2];
[ZZ,IndSort] = sortrows(ZZ);
sizeZZ = size(ZZ,1);
ZTemp = ZZ - [ZZ(sizeZZ,:) ; ZZ(1:(sizeZZ-1),:)];
I = find(sum(abs(ZTemp),2)>0);
INull = find(sum(abs(ZTemp),2)==0);
if isempty(I)                % 12/26/01 -- Identical monomials
    I = 1;
    INull = 2;
end;
Z = ZZ(I,:);

% Constructing permutation matrix
Ind = Ind(IndSort,:);
for i = INull
    Ind(i-1,2) = Ind(i,2);
    Ind(i,2) = sizeZ2+1;
end;
Ind = Ind(I,:);

R1 = [speye(sizeZ1) sparse(sizeZ1,size(I,1)-sizeZ1)];   % 12/12/01
R1 = R1(:,Ind(:,1));
R2 = [speye(sizeZ2) sparse(sizeZ2,size(I,1)-sizeZ2)];    % 12/12/01
R2 = R2(:,Ind(:,2));

Z = makesparse(Z);          % 12/12/01


