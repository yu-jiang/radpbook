function b = subs(a,old,new)
% function B = subs(A,Old,New);
%
% DESCRIPTION
%  Symbolic Substitution.
%
% INPUTS
%   A: Nr-by-Nc polynomial array
%   Old: No-by-1 array of polynomial variables or No-by-1 cell array
%       of characters. The entries of Old must be unique.
%   New: No-by-Npts array of polynomials or doubles.  If Npts>1 then
%       A must be a column or row vector.
%
% OUTPUTS
%   B:  polynomial. B is always returned as a polynomial. Use 'double'
%       to convert B to a double when the final result is a constant.
%       If Npts=1 then B is Nr-by-Nc.  If Npts>1, B is Nr-by-Npts when
%       Nc=1 and Npts-by-Nc otherwise.
%
% SYNTAX
%   B = subs(A,Old,New);
%     Replaces variables in Old with the corresponding entries in New.
%   B = subs(A);
%     Replaces all variables in A with values in the BASE workspace.
%   B = subs(A,New);
%     If New is an 1-by-1 polynomial array then this is equivalent
%     B=subs(A,A.varname{1},New).  Otherwise, this is equivalent to
%     B=subs(A,New(:,1),New(:,2:end)).
%
% EXAMPLE
%  pvar x1 x2 y
%  x=[x1;x2];
%  p=2*(x1+x2)^2+5;
%  subs(p,x,[1;2])
%  subs(p,x,[0 1 1; 1 0 2])
%  subs(p,x1,y)
%
% See also double

% 3/25/2003  PJS Initial Coding
% 8/6/2009   PJS Bug for numeric subs where a is constant (a.nvar==0)
% 12/10/2010 PJS Update syntax. Use collect/peval for partial numeric subs.

%---------------------------------------------------------------------
% Convert to basic syntax: B = subs(A,Old,New);
%---------------------------------------------------------------------
a = polynomial(a);
if nargin == 1
    % Syntax: B = subs(A);
    new = polynomial;
    nva = length(a.varname);
    for i1 = 1:nva
        vari = a.varname{i1};
        if evalin('base', ['exist(''' vari ''')'] )
            val = evalin('base',vari);
            vari = pvar(vari);
            if ~isequal(vari,val)
                new = [new; vari val ];
            end
        end
    end
    if isempty(new)
        b = a;
        return;
    else
        b = subs(a,new);
        return;
    end
elseif nargin == 2
    if size(old,2)>=2
        % Syntax:  B = subs(A,OldNew);
        L.type = '()';
        L.subs = {':',2:size(old,2)};
        new = subsref( old, L);
        L.subs = {':',1};
        old = subsref( old, L);
        b=subs(a,old,new);
        return;
    else
        % Syntax:  B = subs(A,New);
        b=subs(a,a.varname{1},old);
        return;
    end
end

if isa(new,'double') && isvector(new) && all(size(old)==[1 1])
    % Convert column vector to row vector for single var replacement
    new = new(:)';
end
    

%---------------------------------------------------------------------
% Convert to basic dimension case A
%  A) a is Nr x 1, old is No x 1, new is No x Npts --> b is Nr x Npts
%  B) a is 1 x Nc, old is No x 1, new is No x Npts --> b is Npts x Nc
%  C) a is Nr x Nc (Nr,Nc>1), old/new are No x 1 --> b is Nr x Nc
%---------------------------------------------------------------------

% Convert old to a cell array of chars
if ispvar(old)
    old = char(old);
elseif ~iscellstr(old)
    error(['Old must be a No-by-1 array of polynomial variables or '...
        'a No-by-1 cell array of characters.']);
end
if length(old)~=length(unique(old))
    error('The entries of Old must be unique');
end

% Convert constant new to a double
if isa(new,'polynomial')
    if isdouble(new)
        new = double(new);
    end
elseif ~isa(new,'double');
    error('New must be a No-by-Npts array of polynomials or doubles.');
end

old = old(:);
No = length(old);
[Nr,Nc] = size(a);
sznew = size(new);
if sznew(1)~=No
    error('Row dimension of new must equal length of old');
else
    Npts = sznew(2);
end
if Nr>1 && Nc>1 && Npts>1
    error(['Column dimension of new is > 1. A must be a row or '...
        'column vector.'])
end

if Nr>1 && Nc>1
    %  C) a is Nr x Nc (Nr,Nc>1), old/new are No x 1 --> b is Nr x Nc
    L.type = '()'; L.subs = {':'};
    acol = subsref(a,L);
    b = subs( acol, old, new);
    b = reshape(b,[Nr,Nc]);
    return
elseif Nr==1 && Nc>1
    %  B) a is 1 x Nc, old is No x 1, new is No x Npts --> b is Npts x Nc
    b = subs(a',old,new)';
    return
end

%---------------------------------------------------------------------
% Perform subs in basic case:
% a is Nr x 1, old is No x 1, new is No x Npts --> b is Nr x Npts
%---------------------------------------------------------------------

% Get polynomial data
acoef = a.coefficient;
adeg = a.degmat;
avar = a.varname;

% Handle constant a 
if isempty(avar)
    b=repmat(double(a),[1 Npts]);
    return
end

% Perform numeric or symbolic substitution
[ismem,idxa]=ismember(avar,old);
if isa(new,'double') && all(ismem)
    % Numeric replacement of all vars of A-->Use peval for speed.
    % new values passed to peval must be a full matrix
    new = full(new(idxa,:));
    b = peval(new, acoef,adeg);
    b = polynomial(b);
    
elseif isa(new,'double')
    % Numeric replacement of some vars in A
    % Express: p(x,y) = g0(x)+g(x)*h(y) where y are the vars in old
    % Use peval to do the replacement on h(y) and then recombine terms.
    x = avar(~ismem);
    [g0,g,h] = collect(a,x);
    if isempty(h)
        b = repmat(g0,[1 Npts]);
    else
        [ismem,idxh]=ismember(h.varname,old);
        new = full(new(idxh,:));
        hval = peval(new, h.coefficient,h.degmat);
        b = repmat(g0,[1 Npts])+g*hval;
    end
elseif isa(new,'polynomial')
    % Symbolic replacement
    
    % Express: p(x,y) = g0(x)+g(x)*h(y) where y are the vars in old
    x = a.varname(~ismem);
    [g0,g,h] = collect(a,x);
    if isempty(h)
        b = repmat(g0,[1 Npts]);
    else
        % Re-order rows of hdeg to correspond with ordering in h
        hcoef = h.coefficient;
        hdeg = h.degmat;
        [ridx,cidx]=find(hcoef);
        hdeg(cidx,:) = hdeg(ridx,:);

        % Re-order rows of new to align with ordering in h.varname
        [ismem,idxh]=ismember(h.varname,old);
        L.type = '()';
        L.subs = {idxh,':'};
        new = subsref(new,L);

        % Replace vars in h with entries of new
        Nh = size(hdeg,1);
        hval = polynomial( ones(Nh,Npts) );
        for i1=1:size(hdeg,2);
            hi = hdeg(:,i1);
            hi = repmat(hi,[1 Npts]);

            L.subs = {i1, 1:Npts};
            newi = subsref(new,L);
            newi = repmat(newi,[Nh,1]);
            if i1==1
                hval = (newi.^hi);
            else
                hval = hval.*(newi.^hi);
            end
        end

        % Recombine terms
        b = repmat(g0,[1 Npts])+g*hval;
    end
%     % OLD CODE
%     b = polynomial( zeros(Nr,Npts) );
%     L.type = '()';
%     for i1=1:Npts
%         L.subs = {':',i1};
%         newi = subsref(new,L);
%         tmp = LOCAL_symbolicsubs(a,old,newi);
%         b = subsasgn(b,L,tmp);
%     end
end


% %---------------------------------------------------------------------
% % LOCAL_symbolicsubs
% %---------------------------------------------------------------------
% function b = LOCAL_symbolicsubs(a,old,new)
% 
% if 0
%     % XXX ORIGINAL VERSION
%     b=a;
%     ntb = size(b.degmat,1);
%     nvb = size(b.degmat,2);
%     bcoef = b.coefficient;
%     bdeg = b.degmat;
%     bvar = b.varname;
%     szb = size(b);
%     
%     for i1 = 1:length(old);
%         temp = strmatch(old{i1},b.varname,'exact');
%         if isempty(temp)
%             idx(i1) = -1;
%             %error(['Undefined variable ' old{i1}]);
%         else
%             idx(i1) = temp;
%         end
%     end
%     idx2=setdiff(1:nvb,idx);
%     
%     b=polynomial(0);
%     L.type = '()';
%     for i1 = 1:ntb
%         coef = bcoef(i1,:);
%         if isempty(idx2)
%             temppoly = polynomial(reshape(coef,szb));
%         else
%             %coef = coef(:);
%             temppoly = polynomial(coef,bdeg(i1,idx2),bvar(idx2),szb);
%         end
%         
%         for i2 = 1:length(old)
%             if idx(i2)>0
%                 L.subs = {i2};
%                 newi = subsref(new,L);
%                 temppoly = temppoly*power(newi,bdeg(i1,idx(i2)));
%             end
%         end
%         b=b+temppoly;
%     end
% else
%     % XXX PJS 11/12/09 -- Update to limit number of calls to PLUS and
%     % COMBINE based on a moderately large subs example from George Hines.
%     b=a;
%     for i1 = 1:length(old);
%         % Get updated poly data
%         szb = size(b);
%         bcoef = b.coefficient;
%         bdeg = b.degmat;
%         bvar = b.varname;
%         
%         % Find idx for var to be replaced
%         varidx = find( strcmp(old{i1},b.varname) );
%         if ~isempty(varidx)
%             % Get new symbolic expression
%             L.type = '()';
%             L.subs = {i1};
%             newi = subsref(new,L);
%             if ispvar(newi)
%                 % Simply replace old and new var names
%                 b.varname{varidx} = newi.varname{1};
%             else
%                 % Collect terms that don't depend on OLD
%                 degidx = find( bdeg(:,varidx) == 0 );
%                 if isempty(degidx)
%                     b = polynomial(zeros(szb));
%                 else
%                     b = polynomial(bcoef(degidx,:),bdeg(degidx,:),bvar,szb);
%                 end
%                 
%                 % Update terms that depend on OLD
%                 olddeg = unique(bdeg(:,varidx));
%                 olddeg( olddeg == 0 ) = [];
%                 for i2 = 1:length(olddeg)
%                     % Get polynomial associated with newi^deg
%                     temppoly = power(newi,olddeg(i2));
%                     
%                     % Get terms of b that have OLD^deg
%                     degidx = find(bdeg(:,varidx)==olddeg(i2));
%                     coef2 = bcoef(degidx,:);
%                     deg2 = bdeg(degidx,:);
%                     deg2(:,varidx) = [];
%                     var2 = bvar;
%                     var2(varidx) = [];
%                     
%                     % Create poly with OLD^deg replaced by newi^deg
%                     btmp = polynomial(coef2,deg2,var2,szb)*temppoly;
%                     
%                     % Add result to running total
%                     b=b+btmp;
%                 end
%                 
%             end
%         end
%     end
%     
%     
% end


