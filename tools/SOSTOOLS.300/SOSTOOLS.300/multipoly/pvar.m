function varargout = pvar(varargin)
% function p = pvar(varargin)
%
% DESCRIPTION 
%   Create variables (i.e. monomials of degree 1).
%   
% INPUTS 
%   X1,X2,...: Character strings used to name variables.
%
% OUTPUTS  
%   p: pvar
%
% SYNTAX 
%   pvar('x1','x2','x3')  
%   pvar x1 x2 x3  
%     Both of these function calls create monomials of degree 1 in the
%     caller workspace with the given names. Any number of pvars can be 
%     created.
%   p1 = pvar('x1') 
%     Creates a pvars named x1 and assigns it to the output variable p1.
%   [p1,p2,...] = pvar('x1','x2',...) 
%     Creates many pvars and assigns them to the output variables.
%
% See also mpvar

% 6/8/2002: PJS  Initial Coding  
% 12/24/2009: PJS Added output syntax to allow 5*pvar('a')

for i1 = 1:nargin
    if ischar(varargin{i1})
        if nargout==0
            assignin('caller',varargin{i1}, polynomial(varargin(i1)) );
        else
            varargout{i1} = polynomial(varargin(i1));
        end
    else
        error('Inputs must be strings')
    end
end
  