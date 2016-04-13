function B = cleanpoly(A,tol,deg)
% function B = cleanpoly(A,tol,deg)
%
% DESCRIPTION
%   Cleans up the input polynomial.  The output polynomial includes only
%   terms whose coefficients have magnitude greater than or equal to TOL
%   and whose monomial degree is specified by DEG.
%
% INPUTS
%   A: polynomial
%   tol: scalar double specifying the coefficient tolerance
%   deg: vector of non-negative integers specifying the degrees of
%        mononmials to retain. Alternatively deg can be an N-by-2
%        cell array with deg{i,1} specifying a variable and
%        deg{i,2} specifying a vector of non-negative integers.
%        This will retain only monomials whose degree in variable
%        deg{i,1} is specified in deg{i,2}.
%
% OUTPUTS
%   B: polynomial which only contains the terms of A whose coefficients
%      have magnitude greater than or equal to tol and whose monomial
%      degree is listed in deg.
%
% SYNTAX
%   B=cleanpoly(A,tol);
%   B=cleanpoly(A,[],deg);
%   B=cleanpoly(A,tol,deg);
%
% EXAMPLE
%   pvar x1 x2 u;
%   p = 9*u^3 + u*x1^2 + 1e-6*u^2*x1*x2 + 1e-5*u*x2^2 + 2*x1^3 ...
%        - x1*x2 + 3*u + x1 + 2*x2;
%
%   % Remove terms whose coefficients has magnitude < tol
%   tol = 1e-4;
%   p1 = cleanpoly(p,tol)
%
%   % Retain linear and quadratic terms
%   p2 = cleanpoly(p,[],1:2)
%
%   % Retain terms linear in u but of degree 0,1,2,3 in x1 and x2
%   p3 = cleanpoly(p,[],{x1, 0:3; x2, 0:3; u 1})

% 1/28/2008: PJS  Added option to grab terms of a specified degree
% 4/22/2009: PJS  Added functionality to specify deg as a cell array

if nargin<3
    deg = [];
end

if ~isempty(tol)
    if isa(tol,'double') && isscalar(tol) && tol>=0
        Acoef = get(A,'Coefficient');
        idx = find(abs(Acoef) < tol);
        if ~isempty(idx)
            Acoef(idx) = 0;
            chkval = 0; % skip validity check
            A=polynomial(Acoef,A.degmat,A.varname,size(A),chkval);
        end
    else
        error('tol must be a non-negative double');
    end
end

if ~isempty(deg)
    
    if isa(deg,'double') && ndims(deg)==2 && ...
            all(floor(deg)==ceil(deg)) && all(deg>=0)
        
        % deg is a vector of non-negative integers
        
        Adeg = sum(A.degmat,2);
        idx = [];
        for i1=1:length(deg)
            %idx=find(Adeg~=deg(i1));
            idx=[idx; find(Adeg==deg(i1))];
        end
        idx = setdiff(1:length(Adeg),idx);
        
        Acoef = A.coefficient;
        Acoef(idx,:) = 0;
        chkval = 0; % skip validity check
        A=polynomial(Acoef,A.degmat,A.varname,size(A),chkval);
        
    elseif iscell(deg) && ndims(deg)==2 && size(deg,2)==2
        % deg is a cell array of pvars and vectors of non-neg integers
        
        Acoef = A.coefficient;
        for i1=1:size(deg,1)
            var = deg{i1,1};
            if ispvar(var) && all(size(var)==[1 1])
                var = char(var);
            elseif ~ischar(var)
                error('First column of cell array deg must contain strings or pvars');
            end
            vidx = find(strcmp(A.varname,var));
            
            if ~isempty(vidx)
                degi = deg{i1,2};
                if ~( isa(degi,'double') && ndims(degi)==2 && ...
                        all(floor(degi)==ceil(degi)) && all(degi>=0) )
                    error(['Second column of cell array deg must '...
                        'contain a vector of non-negative integers.']);
                end
                
                Adeg = A.degmat(:,vidx);
                cidx = [];
                for i2=1:length(degi)
                    cidx=[cidx; find(Adeg==degi(i2))];
                end
                cidx = setdiff(1:length(Adeg),cidx);
                Acoef(cidx,:) = 0;
            end
        end
        A=polynomial(Acoef,A.degmat,A.varname,size(A));
    else
        error('deg must be a vector of non-negative integers or an Nx2 cell array');
    end
    
end

% Combine to remove terms
B = combine(A);
