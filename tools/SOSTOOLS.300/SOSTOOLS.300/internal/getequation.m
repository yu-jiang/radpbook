function [At,b,Z] = getequation(symexpr,vartable,decvartable,varmat)
%function [At,b,Z] = getequation(symexpr,vartable,decvartable,decvartablename)
%
% GETEQUATION --- Convert a symbolic expression to At, b, and Z
%         used in an SOS program.

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

% 03/01/02 - SP
% 03/10/02 - SP -- Use Maple
% 07/10/13 - JA&GV -- Use the matlab's symbolic engine
% 09/11/13 - PJS -- Addition for multipoly objects

if isa(symexpr,'polynomial')
    % PJS: Handle Polynomial Objects including the Matrix Case
    decvarnum = numel(decvartable);
    vartable = [vartable; varmat];
    cvartable = char(vartable);
    dimp = size(symexpr,2);
    
    % Collect symexpr = g(c)*h(x) where x are the vars in vartable,
    % c are the decision variables, and h(x) is a vector of monoms.
    % in the actual polynomial vars.  Use this to find unique
    % monomials in the polynomial variables.
    [g0,g,h] = collect(symexpr(:),setdiff(symexpr.varname,cvartable));
    g = g(:);
    if ~isequal(g0,0)
        g = [g0;g];
        h = [1; h];
    end
    [nmon,nvar] = size(h.degmat);
    
    % Reorder the monomials matrix with variables in the order
    % listed in cvartable and sorted monomials
    Z = zeros(nmon,nvar);
    [~,idx]=ismember(h.varname,cvartable);
    Z(:,idx) = full(h.degmat);
    Z = sortNoRepeat([],Z);
    
    % Initialize coefficient variables
    [nmon,nvar] = size(Z);
    coeffnts = sparse(nmon*dimp,dimp);
    if decvarnum>0
        coeffnts_decvar = polynomial(sparse(nmon*dimp,dimp));
    end
    
    % Define coefficients
    for i = 1:dimp
        for j = i:dimp
            exprij = symexpr(i,j);
            [g0,g,h] = collect(exprij,setdiff(exprij.varname,cvartable));
            g = g(:);
            if ~isequal(g0,0)
                g = [g0;g];
                h = [1; h];
            end
            coefmatr = double(subs(g,decvartable,zeros(decvarnum,1)));
            monmatr = zeros(length(h),nvar);
            [~,idx]=ismember(h.varname,cvartable);
            monmatr(:,idx) = h.coef'*h.degmat;
            
            for k = 1:length(h);
                s_ijk = coefmatr(k,1);
                s_ijk_decvar = g(k);
                mon_k= monmatr(k,:);
                
                [val,ind_k] = max(sum((Z == kron(ones(nmon,1),mon_k)),2));
                coeffnts((ind_k-1)*dimp+i,j) = s_ijk;
                coeffnts((ind_k-1)*dimp+j,i) = s_ijk;
                
                coeffnts_decvar((ind_k-1)*dimp+i,j) = s_ijk_decvar;
                coeffnts_decvar((ind_k-1)*dimp+j,i) = s_ijk_decvar;
            end 
        end
    end
    
    % Constructing At and b matrices
    At = sparse(decvarnum,size(Z,1)*dimp^2);    
    if decvarnum>0
        for i=1:nmon
            Mivec = reshape(coeffnts_decvar((i-1)*dimp+1:i*dimp,:),dimp^2,1);
            At(:,(i-1)*dimp^2+1:i*dimp^2) = sparse(-double(jacobian(Mivec,decvartable))');
        end        
    end
    b = sparse(coeffnts);
    
else
    % PJS: Original Code to Handle Symbolic Objects
    
    if strcmp(decvartable,'[]');
        decvarnum = 0;
    else
        decvarnum = length(find(decvartable == ','))+1;
    end;
    expr = evalin(symengine,symexpr);
    %vartable = evalin(symengine,vartable);
    if  length(varmat)==2
        vartable = evalin(symengine,vartable);
    else
        vartable = evalin(symengine,[vartable(1:end-1),',',varmat(2:end)]);  %JA&GV 07/06/13
    end
    
    decvartable = evalin(symengine,decvartable);
    charvartable = converttochar([vartable]);
    
    FPexpr = feval(symengine,'expand',expr);
    
    %JA&GV 07/13 commented the below to implement the matrix case
    %FPexpr = feval(symengine,'collect',FPexpr,charvartable);
    
    dimp = size(FPexpr,2);%JA&GV jul/04/2013 sets the dimension of the matrix of polynomials
    
    for i = 1:size(FPexpr,1)
        for j = 1:dimp
            FPexpr(i,j) = feval(symengine,'collect',FPexpr(i,j),charvartable);
        end
    end
    
    
    
    if isempty(decvartable)==0%JA&GV 07/13 the code below  addresses the matrix case
        
        
        Zfull = [];
        for i = 1:dimp
            for j = i:dimp
                coefmon = feval(symengine,'poly2list',FPexpr(i,j),charvartable);
                coefmonmatr = subs(coefmon.',decvartable,zeros(1,length(decvartable)));
                Z = double(coefmonmatr(:,2:end));
                Zfull = sortNoRepeat(Zfull,Z);
            end
        end
        Z = Zfull;
        [nmon,nvar] = size(Z);
        coeffnts = sparse(nmon*dimp,dimp);
        coeffnts_decvar = sym(sparse(nmon*dimp,dimp));
        
        
        for i = 1:dimp
            for j = i:dimp
                coefmon = feval(symengine,'poly2list',FPexpr(i,j),charvartable);
                coefmonmatr = subs(coefmon.',decvartable,zeros(1,length(decvartable)));
                for k = 1:size(coefmonmatr,1)
                    s_ijk = coefmonmatr(k,1);
                    
                    dummy_decvar = reshape(coefmon(k),2,1);
                    s_ijk_decvar = dummy_decvar;
                    
                    mon_k= double(coefmonmatr(k,2:end));
                    [val,ind_k] = max(sum((Zfull == kron(ones(nmon,1),mon_k)),2));
                    coeffnts((ind_k-1)*dimp+i,j) = s_ijk;
                    coeffnts((ind_k-1)*dimp+j,i) = s_ijk;
                    
                    coeffnts_decvar((ind_k-1)*dimp+i,j) = s_ijk_decvar;
                    coeffnts_decvar((ind_k-1)*dimp+j,i) = s_ijk_decvar;
                    
                    
                end
                
            end
        end
        
        
    else%JA&GV 07/2013 the code below  addresses the matrix case
        
        
        Zfull = [];
        for i = 1:dimp
            for j = i:dimp
                coefmon = feval(symengine,'poly2list',FPexpr(i,j),charvartable);
                coefmon = coefmon.';
                for k = 1:length(coefmon)
                    dummyvar = reshape(coefmon(k),2,1);
                    Z(k,:) = double(dummyvar(2));
                end
                Zfull = sortNoRepeat(Zfull,Z);
                Z = [];
            end
        end
        Z = sparse(Zfull);
        [nmon,nvar] = size(Z);
        coeffnts = sparse(nmon*dimp,dimp);
        
        for i = 1:dimp
            for j = i:dimp
                coefmon = feval(symengine,'poly2list',FPexpr(i,j),charvartable);
                coefmon = coefmon.';
                for k = 1:length(coefmon)
                    dummyvar = reshape(coefmon(k),2,1);
                    s_ijk = double(dummyvar(1));
                    mon_k= double(dummyvar(2));
                    [val,ind_k] = max(sum((Zfull == kron(ones(nmon,1),mon_k))'));
                    coeffnts((ind_k-1)*dimp+i,j) = s_ijk;
                    coeffnts((ind_k-1)*dimp+j,i) = s_ijk;
                end
            end
        end
        
    end
    
    % Constructing At and b matrices
    At = sparse(decvarnum,size(Z,1)*dimp^2);
    
    %JA&GV may/2013 uses the jacobian function to derive the expression
    if isempty(decvartable)==0
        for i = 1:size(Z,1)
            Mivec = reshape(coeffnts_decvar((i-1)*dimp+1:i*dimp,:),dimp^2,1);
            At(:,(i-1)*dimp^2+1:i*dimp^2) = sparse(-double(jacobian(Mivec,decvartable))');
        end
    end
    
    b = sparse(coeffnts);
end


