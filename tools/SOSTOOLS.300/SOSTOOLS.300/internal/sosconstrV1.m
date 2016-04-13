function sos = sosconstr(sos,Type,symexpr)
% SOSCONSTR --- Add a new constraint (equality/inequality)
%    to an SOS program
%
% SOSP = sosconstr(SOSP,TYPE,EXPR)
%
% SOSP is the sum of squares program.
% The new constraint is described by TYPE as follows:
%   TYPE = 'eq'   : constraint is an equality, viz., f(x) = 0
%   TYPE = 'ineq' : constraint is an inequality, viz., f(x) >= 0 (an SOS)
% EXPR is the expression in the left hand side of the constraint, i.e., f(x).
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


% 12/24/01 - SP
% 03/01/02 - SP -- New syntax
% 04/06/03 - PJS -- Handle poly objects w/out for-loops

sos.expr.num = sos.expr.num+1;
expr = sos.expr.num;
sos.expr.type{expr} = Type;


if isfield(sos,'symvartable')
    charvartable = converttochar([sos.vartable]);
    isvectorvar = size(symexpr,2)==1;
    for i = 1:size(symexpr,1)
        if isvectorvar
            degcheck = evalin(symengine,char(symexpr(i)));
        else
            degcheck = evalin(symengine,char(symexpr(i,i)));
        end
        degcheck = feval(symengine,'expand',degcheck);
        degcheck = feval(symengine,'collect',degcheck,charvartable);
        degcheckmon = feval(symengine,'poly2list',degcheck,charvartable);
        mondeg = zeros(length(degcheckmon),1);
        for j = 1:length(degcheckmon)
            test = degcheckmon(j);
            mondeg(j) = sum(double(test(2)));
        end
        mondegmax = max(mondeg);
        mondegmin = min(mondeg);
        if mod(mondegmax,2)~=0||mod(mondegmin,2)~=0
            error('SOSTOOLS: Leading and trailing degree of SOS expression can not be odd. I suggest you check your polynomial expression.');
        end
    end
    
    
    [sos.expr.At{expr},sos.expr.b{expr},sos.expr.Z{expr}] = ...
        getequation(char(symexpr),sos.vartable,sos.decvartable,sos.varmat.vartable);
else
    
    if 0 %isscalar(symexpr)
        % PJS 9/9/2013: Original Code for scalar SOS Constraint
        
        % Pull out information from the polynomial
        coef = symexpr.coefficient;
        deg = symexpr.degmat;
        var = symexpr.varname;
        
        % Sort variables: decision stacked on non-decision
        [dummy,idxdvar1,idxdvar2] = intersect(var,sos.decvartable);
        ldv = length(idxdvar1);
        [dummy,idxvar1,idxvar2] = intersect(var,sos.vartable);
        var = var([idxdvar1; idxvar1]);
        deg = deg(:,[idxdvar1; idxvar1]);
        
        % Sort terms: Stack terms without decision variables above
        % terms with decision variables.
        [deg,sortidx] = sortrows(deg);
        coef = coef(sortidx);
        if ldv~=0
            [idx,dummy1]=find( deg(:,1:ldv) );
            nvterm = min(idx)-1;
        else
            nvterm = size(deg,1);
        end
        D12 = deg(1:nvterm,ldv+1:end);
        D22 = deg(nvterm+1:end,ldv+1:end);
        D21 = deg(nvterm+1:end,1:ldv);
        
        % Get the equality constraints
        if any( sum(D21,2) > 1 )
            error(['The expression is not linear in the decision' ...
                ' variables']);
        else
            
            % Extract monomials
            [Ztemp,idx1,idx2] = unique([D12; D22],'rows');
            lZ = size(Ztemp,1);
            sos.expr.Z{expr} = sparse(lZ,length(sos.vartable));
            sos.expr.Z{expr}(:,idxvar2) = flipud(Ztemp);
            
            % Form the equality constraints
            Ac = sparse(length(sos.decvartable),lZ);
            b = sparse(lZ,1);
            
            % If symexpr has not been combined down, need to check that
            % idx2(1:nvterm) are all unique....
            b(idx2(1:nvterm)) = coef(1:nvterm);
            if ldv~=0
                [dv,dummy]=find(D21');
                tempidx = sub2ind(size(Ac), idxdvar2(dv) , idx2(nvterm+1:end) );
                Ac(tempidx) = coef(nvterm+1:end);
            end
            
            sos.expr.At{expr} = -fliplr(Ac);
            sos.expr.b{expr} = flipud(b);
            
            
        end
    else
        % PJS 9/9/2013: New Code for matrix SOS inequality
        % Note: Mimic symbolic code in getequation
        decvartable = sos.decvartable;
        if isempty(decvartable)
            decvarnum = 0;
        else
            decvarnum = numel(decvartable);
        end;
        
        % Pull out information from the polynomial        
        coef = symexpr.coefficient;
        deg = symexpr.degmat;
        var = symexpr.varname;
        dimp = size(symexpr,1);

        % Remove decision variables
        [~,idx]=intersect(var,decvartable)
        var(idx) = [];
        deg(:,idx) = [];
        
        Zfull = [];
        for i = 1:dimp
            for j = i:dimp
                idx = sub2ind([dimp dimp],i,j);
                cij = coef(:,idx);
                Zfull = unique([Zfull; deg(find(cij),:)], 'rows');
            end
        end
        Zfull = flipud(Zfull);
        Z = Zfull;
        [nmon,nvar] = size(Z);
        
        coeffnts = sparse(nmon*dimp,dimp);
        if ~isempty(decvartable)
            coeffnts_decvar = polynomial(zeros(nmon*dimp,dimp));
        end
        
        for i = 1:dimp
            for j = i:dimp
                idx = sub2ind([dimp dimp],i,j);
                cij = coef(:,idx);
                nzidx = find(cij);
                
                for k = 1:numel(nzidx)
                    s_ijk = cij(nzidx(k));
                    mon_k = deg(nzidx(k),:);                    
                    [val,ind_k] = max(sum((Zfull == kron(ones(nmon,1),mon_k)),2));
                    coeffnts((ind_k-1)*dimp+i,j) = s_ijk;
                    coeffnts((ind_k-1)*dimp+j,i) = s_ijk;
                    
                    if 1 %~isempty(decvartable)
                        % XXX TO FIX
                        dummy_decvar = reshape(coefmon(k),2,1);
                        s_ijk_decvar = dummy_decvar;
                        coeffnts_decvar((ind_k-1)*dimp+i,j) = s_ijk_decvar;
                        coeffnts_decvar((ind_k-1)*dimp+j,i) = s_ijk_decvar;
                    end
                end
            end
        end
        
        % Constructing At and b matrices
        At = sparse(decvarnum,size(Z,1)*dimp^2);
        
        %JA&GV may/2013 uses the jacobian function to derive the expression
        if ~isempty(decvartable)
            % XXX TO FIX
            for i = 1:size(Z,1)
                Mivec = reshape(coeffnts_decvar((i-1)*dimp+1:i*dimp,:),dimp^2,1);
                At(:,(i-1)*dimp^2+1:i*dimp^2) = sparse(-double(jacobian(Mivec,decvartable))');
            end
        end
        b = sparse(coeffnts);
        
        % Store constraints
        sos.expr.At{expr} = At;
        sos.expr.b{expr} = b;
        sos.expr.Z{expr} = Z;
    end
    
end



