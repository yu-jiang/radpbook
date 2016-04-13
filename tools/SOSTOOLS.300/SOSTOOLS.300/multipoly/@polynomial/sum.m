function varargout=sum(A,dim)
% function C=sum(A,Dim)
%
% DESCRIPTION
%   Sum the elements of a polynomial vector or matrix.
%
% INPUTS
%   A: polynomial vector or matrix
%   Dim: dimension along which summation is to be performed
%      (Default: Dim = 1)
%
% OUTPUTS
%   C: polynomial or polynomial vector, the result of summation
%
% SYNTAX
%   C = sum(A)
%     If A is a vector then C is the sum of the elements in A.  If A is a
%     matrix then C is a row vector with the sum over each column of A.
%   C = sum(A,1)
%     C is a row vector with the sum over each column of A.
%   C = sum(A,2)
%     C is a column vector with the sum over each row of A.

% Initial coding: 01/31/03 -- SP
% 3/25/03 PJS   If x is a column vector, sum(x,2) should return x.
%               If nargout = 0, display the result.
% 11/12/09 PJS  Removed for loops in vector case for speed. Have matrix
%               case call vector case in a loop
% 11/23/10 PJS  Vectorized matrix case

% Promote A
A = polynomial(A);
szA = size(A);
if (any(szA == 1) && nargin==1) || ( szA(2)==1 && dim==1 ) ...
        || ( szA(1)==1 && dim==2 )
    % Vectorized for speed
    C = A;
    C.coefficient = sum(C.coefficient,2);
    C.matdim = [1 1];
    C = combine(C);
elseif nargin == 1 || dim == 1
    % Call dim=2 vectorized code
    C = sum(A',2)';
    
%     C = polynomial( zeros(1,szA(2)) );
%     for i1 = 1:size(A,2)
%         s.type = '()';
%         s.subs = {':',i1};
%         Ai = subsref(A,s);
%         
%         s.subs = {i1};
%         C =  subsasgn(C,s,sum(Ai));
%     end
else
    % Vectorized for speed
    C = A;    
    nta = size(C.coefficient,1);

    tmp = reshape(C.coefficient,[nta*szA(1) szA(2)]);
    tmp = sum( tmp , 2 );    
    C.coefficient = reshape( tmp, [nta, szA(1)]);
    C.matdim = [szA(1) 1];
    C = combine(C);
    
%     C = polynomial( zeros(szA(1),1) );
%     for i1 = 1:size(A,1)
%         s.type = '()';
%         s.subs = {i1,':'};
%         Ai = subsref(A,s);
%         
%         s.subs = {i1};
%         C =  subsasgn(C,s,sum(Ai));
%     end
end

varargout{1} = C;

