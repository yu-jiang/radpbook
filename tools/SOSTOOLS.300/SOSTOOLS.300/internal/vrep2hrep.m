function [A,B] = vrep2hrep(pts);
% 
% Given a set of points, computes an H-representation of its convex hull  
%
%  A x <= B
%
% Calls Komei Fukuda's CDD (floating point version)

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

% PP 22/10/2002  
  
[npts,dim]=size(pts);  
  
cdd = cddpath;

filename = 'temp.ext' ;

fid = fopen(filename,'w');

fprintf(fid,'*\n');
fprintf(fid,'* Generated automatically by vrep2hrep.m\n');
fprintf(fid,'*\n');
fprintf(fid,'V-representation\n');
fprintf(fid,'begin\n');
fprintf(fid,'%d %d real\n',npts,dim+1);

fmt = ['1 ' repmat('%.15E ',1,dim) '\n'];

% Print the data
fprintf(fid,fmt,pts');

fprintf(fid,'end\n');
fprintf(fid,'stdout_off\n');
fclose(fid);

% Now, run CDD
cmd = [cdd ' ' filename];

vv = version;
if str2num(vv(1))>5
  [dummy1,dummy2] = system(cmd);
else
  if ~isunix
    [dummy1,dummy2] = dos(cmd);
  else
    [dummy1,dummy2] = unix(cmd);
  end;  
end;

% Read the output

filename = 'temp.ine' ;

fid = fopen(filename,'r');

% Skip everything, until 'begin'
tline = [];
while strcmp(tline,'begin') == 0;
 tline = fgetl(fid) ;
end

nineqs = fscanf(fid,'%d',1);
dims   = fscanf(fid,'%d',1);
dummy2 = fscanf(fid,'%s',1);     

H = fscanf(fid,'%E',[dims,nineqs]);
fclose(fid);

% Get rid of the files
delete('temp.ext');
delete('temp.ine');

H = H' ;

B = H(:,1);
A = -H(:,2:dims);
