function c = blkdiag(varargin)
% function C=blkdiag(A,B);
%
% DESCRIPTION
%   Block diagonal concatenation of polynomial matrices.
%
% INPUTS
%   A,B: polynomials
%
% OUTPUTS
%   C:  block diagonal concatenation of input matrices.
%
% SYNTAX
%   C = blkdiag(A,B,...);                       [A 0 ...]
%       Block diagonal concatenation produces   [0 B ...]
%                                               [0 0 ...]
% See also diag, horzcat, vertcat

% 6/8/2002: PJS  Initial Coding

if nargin==1
    c = varargin{1};
else
    % Promote a to polynomial
    a = polynomial(varargin{1});
    sza = size(a);
    
    % Promote b to polynomial
    b = polynomial(varargin{2});
    szb = size(b);
    
    if isempty(b);
        c = a;
    elseif isempty(a);
        c = b;
    else
        t12 = polynomial( zeros(sza(1),szb(2)) );
        t21 = polynomial( zeros(szb(1),sza(2)) );
        c = [a t12; t21 b];
    end
    
    if nargin>2
        c = blkdiag(c,varargin{3:end});
    end
end

