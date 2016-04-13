function b = sprepmat(a,m,n)
% function B = sprepmat(A,M,N)
%
% DESCRIPTION 
%   B is a sparse matrix consisting of an MxN tiling of the 
%   sparse matrix A.
%   [repmat fails on scalar, sparse matrices in Matlab 6.5:
%     >>repmat(sparse(5),2,3) 
%    This yields an error.]
%   
% INPUTS 
%   A: sparse matrix
%   M,N: Copies of A in the row and column directions
%
% OUTPUTS  
%   B: sparse matrix
%  
% SYNTAX 
%   B = sprepmat(A,M,N);

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

% 1/30/2003: PJS  Initial Coding  

[nra,nca] = size(a);
ridx = (1:nra)';
ridx = ridx(:,ones(1,m));
cidx = (1:nca)';
cidx = cidx(:,ones(1,n));
b = a(ridx,cidx);
