function [sos] = sosmatrixineq(sos,fM,option)
% SOSMATRIXINEQ --- Creates a SOS constraint from a matrix inequality
% constraint
%
% [SOSP] = sosmatrixineq(SOSP,fM)
%
% SOSP is the sum of squares program.
% fM is a polynomial matrix used to generate the polynomial inequality y'fMy>0

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


% AP - 16/4/2013
% JA - 6/6/2013

if nargin<3
    option='quadraticMineq';%sets the default: asks for sos expression v'M(x)v
end

[n,m] = size(fM);

if n~=m
    disp('ERROR: Matrix fM in inequality fM>0 must be square.');
    return
end


if isfield(sos,'symvartable')
    % Original Code

    if strcmp(option,'quadraticMineq')

        %creates the vector of variables Mvar to generate the quadratic expression
        %M_var'*fM*Mvar
        if n>sos.varmat.count
            s1 = 'syms ';
            for i=sos.varmat.count+1:n
                s1 = strcat(s1,sprintf(' Mvar_%d',i));
            end
            %eval(expression) evaluates the MATLAB code in the string expression.
            eval(s1);
            
            
            Mvar = sym(zeros(n-sos.varmat.count,1));
            for i=sos.varmat.count+1:n
                Mvar(i) = eval(['Mvar_',int2str(i)]);
            end
            
            %updates the vartable in the sos program
            sos.varmat.symvartable = [sos.varmat.symvartable; Mvar];
            Mvarctable = sym2chartable(Mvar);
            if sos.varmat.count==0
                sos.varmat.vartable = [Mvarctable(1:end)];
            else
                sos.varmat.vartable = [sos.varmat.vartable(1:end-1),',',Mvarctable(2:end)];
            end
            sos.varmat.count = n;
        end
        
        varMconst = sos.varmat.symvartable(1:n);
        
        %create the sosconstraint using the sparse multipartite option since it is
        %homogeneous in varMconst
        sos = sosineq(sos,varMconst.'*fM*varMconst,'sparsemultipartite',{sos.symvartable,varMconst}); %GV&JA 6/12/2013
    
    elseif strcmp(option,'Mineq')
        sos = sosineq(sos,fM); %GV&JA 10/01/2013
    end
else
    if strcmp(option,'quadraticMineq')
    % Multipoly Code: PJS 9/9/2013
    if n>sos.varmat.count
        Mvar = polynomial(zeros(n-sos.varmat.count,1));
        for i=sos.varmat.count+1:n
            Mvar(i) = pvar(['Mvar_' int2str(i)]);
        end
              
        %updates the vartable in the sos program
        sos.varmat.vartable = [sos.varmat.vartable; Mvar];
        sos.varmat.count = n;
    end    
    
    %create the sosconstraint using the sparse multipartite option since it is
    %homogeneous in varMconst
    sos = sosineq(sos,Mvar.'*fM*Mvar,'sparsemultipartite',{sos.vartable,Mvar}); % PJS 9/9/2013
    
    elseif strcmp(option,'Mineq')
        sos = sosineq(sos,fM); %GV&JA 10/01/2013
    end
    
end
