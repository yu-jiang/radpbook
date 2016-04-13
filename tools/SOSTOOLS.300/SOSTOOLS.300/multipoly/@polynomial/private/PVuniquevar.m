function [newdegmat,uniquevar] = PVuniquevar(degmat,varname)
% function [NewDegmat,UniqueVarname] = PVuniquevar(Degmat,Varname)
%
% DESCRIPTION
%   (INTERNAL FUNCTION)
%   Find the unique set of variables in varname and form
%   degmat in terms of this unique set.
%
% INPUTS
%   Degmat: degree matrix for a polynomial
%   Varname: variable list for a polynomial
%
% OUTPUTS
%   NewDegmat: degree matrix in terms of unique variable list
%   UniqueVarname: list of unique variables
%
% SYNTAX
%   [NewDegmat,UniqueVarname] = PVuniquevar(Degmat,Varname)

% 11/25/2002: PJS  Initial Coding
%  1/30/2002: PJS  Combine columns using sparse multiply

if 1
    % Find repeated variables:
    %   sortidx: index such that varsort = varname(sortidx)
    %   repeatidx: nonzero indices indicate repeats
    %  rvec: map entries of varsort to numbers 1,..,nuv
    [varsort,sortidx] = sortrows(char(varname));
    repeatidx = all(varsort(2:end,:) == varsort(1:end-1,:),2);
    rvec = cumsum([1; ~repeatidx]);
    
    % Create unique set of variables
    uniqueidx = sortidx;
    uniqueidx( repeatidx ) = [];
    uniquevar = varname(uniqueidx);
    
    % Form degree matrix in terms of new variable set
    %  (Use sparse multiply, degmat*summat, to sum the rows)
    nv = size(degmat,2);
    nuv = length(uniquevar);
    summat = spalloc(nv,nuv,nv);
    idx = sub2ind([nv nuv],sortidx,rvec);
    summat(idx) = 1;
    newdegmat = degmat*summat;
else
    % XXX - This is faster but the ordering of the monomials will depend
    % on how the polynomial was constructed because the ordering of
    % varnames won't be alphabetical when we do the degmat sort.
    
    % Try to retain ordering in vars as much as possible.
    % In many cases this will result in faster sorting of degmat in
    % PVuniqueterm because the blocks of degmat may have already been
    % sorted.  varnames can be alphabetically sorted at the end of
    % combine (if desired).
    
    % Find repeated variables:
    %   sortidx: index such that varsort = varname(sortidx)
    %   repeatidx: indices of repeats in sorted list
    %   delidx: indices of repeats in unsorted list
    [varsort,sortidx] = sortrows(char(varname));
    tmp = [0; all(varsort(2:end,:) == varsort(1:end-1,:),2)];
    repeatidx = find(tmp);
    delidx = sortidx(repeatidx);
    
    % Create unique set of vars (trying to maintain original ordering)
    uniquevar = varname;
    uniquevar( delidx  ) =[];
    
    % Find mapping from varname to uniquevar
    %   tmpidx: vector of indices that relates the repeats in varsort to a
    %     unique set of vars in varsort.  tmpidx is such that variable
    %     i in varsort is the same as variable tmpidx(i) in varsort
    %   v2uidx: mapping from varname to uniquevar
    nv = length(varname);
    tmpidx=zeros(nv,1);
    tmpidx(1) = 1;
    for i1=2:nv
        tmpidx(i1) =  (tmp(i1)==0)*i1 + (tmp(i1)==1)*tmpidx(i1-1);
    end
    v2uidx(sortidx) = sortidx(tmpidx);
    
    % Form degree matrix in terms of new variable set
    %  (Use sparse multiply, degmat*summat, to sum the rows)
    nuv = length(uniquevar);
    summat = spalloc(nv,nuv,nv);
    idx = sub2ind([nv nuv],1:nv,v2uidx);
    summat(idx) = 1;
    newdegmat = degmat*summat;
end



