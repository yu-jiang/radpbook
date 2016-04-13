function varargout = size(a,dim)
% function out=size(A,Dim)
%
% DESCRIPTION 
%   Size of a polynomial matrix.
%   
% INPUTS 
%   A: polynomial 
%   Dim: dimension
%
% OUTPUTS  
%   Out: Size of the polynomial matrix or length
%      along the specified dimension.
%  
% SYNTAX 
%   Out=size(A);
%     Returns the size (1x2 vector) of the polynomial matrix.
%   Out=size(A,Dim);
%     Returns the row dimension (Dim=1) or the column dimension
%     (Dim=2) of the polynomial matrix.
%   [Out1,Out2]=size(A);  
%     Returns the row (Out1) and column (Out2) dimensions as 
%     separate outputs.
  
% 6/8/2002: PJS   Initial Coding
% 6/10/2002: PJS  Size of empty polys = 0
% 11/5/2002: PJS  Can call 'size' with 2 outputs

% Error Checking
  
% Compute Size  
if nargin ==1      
    out = a.matdim;
    if nargout==0 || nargout==1
        varargout{1}=out;
    elseif nargout==2
        varargout{1} = out(1);
        varargout{2} = out(2);
    end  
elseif nargin == 2
    if ~( isa(dim,'double') && isreal(dim) &&  all(size(dim)==[1 1]) &&...
            all(ceil(dim)==floor(dim)) )
        error('Dimension must be a positive integer scalar');
    else
        if dim==1 || dim==2
            out = a.matdim(dim);
        elseif dim>0
            out = 1;
        else
            error('Dimension must be a positive integer scalar');
        end
    end
    varargout{1} = out;  
end  



