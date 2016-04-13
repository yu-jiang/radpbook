function out = set(a,varargin)
% function Out = set(A,Property,Value)
%
% DESCRIPTION
%   Set polynomial properties.
%
% INPUTS
%   A: polynomial
%   Property: polynomial property
%   Value: value of the property
%
% OUTPUTS
%   Out: new polynomial
%
% SYNTAX
%   Out = set(A,Property1,Value1,Property2,Value2,...);

% 6/7/2002: Written by Packard

GSa = PVget(a,'PropNames');
if nargin>=2
    if floor((nargin-1)/2)==ceil((nargin-1)/2)
        for i=1:(nargin-1)/2
            [idx,fstr] = PVcipmatch(GSa.SPropNames,varargin{2*i-1});
            if ~isempty(idx)
                a = PVset(a,GSa.SPropNames{idx},varargin{2*i});
            else
                error(['Property ' varargin{2*i-1} ' is not settable']);
                %warning(['Property ' varargin{2*i-1} ' not found']);
            end
        end
        [aflag,errormsg] = isvalid(a);
        if aflag==1
            a = combine(a);
            if nargout==1
                out = a;
            else
                assignin('caller',inputname(1),a);
            end
        else
            error(errormsg);
        end
    else
        error('Invalid Property/Value list');
    end
else
    disp(GSa.SPropNames);
end



