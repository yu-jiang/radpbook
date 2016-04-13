function varargout = mpvar(var,n,m,opt)
% function P = mpvar(cstr,N,M,opt);
%
% DESCRIPTION
%   Create a polynomial matrix or vector variable
%
% INPUTS
%   cstr: Character string to be used in creating the coefficient vector.
%   N,M: row and column dimensions of polynomial matrix.
%   opt: If N==M, then set opt = 's' to generate a symmetric
%        matrix variable. 
%
% OUTPUTS
%   P: polynomial matrix
%
% SYNTAX
%   P = mpvar('c',N)
%       Creates an NxN polynomial matrix with entries c_i_j.
%   P = mpvar('c',N,M)
%   P = mpvar('c',[N,M])
%       Creates an NxM polynomial matrix with entries c_i_j.
%   P = mpvar('c',N,N,'s')
%       Creates an NxN symmetric polynomial matrix with entries c_i_j.
%   P = mpvar('c',[N,1])
%   P = mpvar('c',[1,N])
%       Creates an Nx1 or 1xN polynomial vector with entries c_i if 
%       N>1.  If N=1 then this creates a pvar named c.
%   mpvar(cstr,N,M)
%      Equivalent to calling eval([cstr '=mpvar(cstr,N,M);']).
%
% EXAMPLE
%   P = mpvar('p',[2,3])
%
% See also pvar

% 11/10/2002 PP    Initial Coding
% 11/10/2002 PJS   Minor modifications for speed and allow
%                  calling sequence with no outputs
% 12/09/2009 PJS   Modified symmetric call for speed
% 11/07/2010 PJS   Added syntax for [N,M] and for vectors

% Argument checking
error(nargchk(2,4,nargin))
error(nargchk(0,1,nargout))
if nargin==2
    m = [];
    opt = [];
elseif nargin==3
    if ischar(m)
        opt = m;
        m = [];
    else
        opt = [];
    end
end
if isempty(m)
    if isscalar(n)
        m = n;
    else
        m = n(2);
        n = n(1);
    end
end

if n~=m && strcmp(opt,'s')
    error('Symmetric option only valid for square matrices');
end

% Coefficient matrix
nt = n*m;
coefficient = speye(nt);

% Degree matrix
degmat = speye(nt);

% Variable names
stridx = int2str( (1:max([n m]))' );
stridx = strjust( stridx ,'left');
stridx = cellstr(stridx);
varname = cell(nt,1);
if n==1 && m==1
    % Scalar
    varname{1} = var;
elseif n>1 && m>1
    % Matrix
    for i1=1:n ;
        tempvar = [var '_' stridx{i1} '_'];
        for i2 =1:m;
            varname{i1+(i2-1)*n} = [tempvar stridx{i2}];
        end
    end
else
    % Vector
    for i1=1:nt
        varname{i1} = [var '_' stridx{i1}];
    end    
end

% Create matrix dimension
matdim = [n m];
chkval = 0; % skip validity check
P = polynomial(coefficient,degmat,varname,matdim,chkval);
if strcmp(opt,'s')
    % Single indices into lower/upper parts of a symmetric matrix
    % sorted so that single index for (i,j) is aligned with (j,i)
    M = reshape(1:n^2,[n n]);
    Ml = tril(M,-1);
    Ms = Ml+Ml';
    [junk,idx]=sort(Ms(:));
    idx=reshape(idx(n+1:end),[2,n*(n-1)/2]);
    lidx = idx(1,:);
    uidx = idx(2,:);
        
    %Create symmetric matrix of pvars
    % [Working with degmat/varname is faster than indexing into P]
    d2 = degmat;
    d2(:,uidx) = d2(:,lidx)+d2(:,uidx);
    d2(:,lidx) = [];
    v2 = varname;
    v2(lidx) = [];
    P = polynomial(coefficient,d2,v2,matdim,chkval);
    P = combine(P);
end

if nargout==0
    assignin('caller', var, P);
else
    varargout{1} = P;
end
