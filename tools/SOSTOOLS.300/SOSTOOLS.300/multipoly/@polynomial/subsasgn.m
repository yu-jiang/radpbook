function b = subsasgn(a,L,RHS)
% function B =  subsasgn(A,L,RHS)
%
% DESCRIPTION
%   Subsassign for polynomial objects.
%
% INPUTS
%   A: polynomial
%   L: a structure array with the fields:
%    type -- string containing '()'  or '.' specifying the
%              subscript type.
%    subs -- Cell array or string containing the actual subscripts.
%   RHS: Value to be assigned.
%
% OUTPUTS
%   B: object after subsassignment
%
% SYNTAX
%   B =  subsasgn(A,L,RHS)

% 6/9/2002: PJS  Initial Coding
% 1/22/2008: PJS Fixed bug to handle a(2,:) =x1 for a undefined or 0-by-0

a = polynomial(a);
sza = size(a);
switch L(1).type
    
    case '.'
        if length(L) == 1
            temp = RHS;
        else
            temp = subsref(a,L(1));
            temp = subsasgn(temp,L(2:end),RHS);
        end
        b = set(a,L(1).subs,temp);
        
    case '()'
        
        % Peform all subsasgn but L(1)
        if length(L)==1
            temp = polynomial(RHS);
        else
            temp = subsref(a,L(1));
            temp = subsasgn(temp,L(2:end),RHS);
        end
        
        %  Three '()'-subsasgn cases
        if length(L(1).subs)==1 && strcmp(L(1).subs{1},':')
            b = PVsubsasgn_colon(a,L,temp);
        elseif length(L(1).subs)==1
            b = PVsubsasgn_1idx(a,L,temp);
        else
            if strcmp(L(1).subs{1},':')
                % To handle the case a(:,1)=x1 when a is undefined or 0-by-0
                if all(sza==[0 0])
                    sza(1) = 1;
                end
                L(1).subs{1} = 1:sza(1);
            end
            if strcmp(L(1).subs{2},':')
                % To handle the case a(:,1)=x1 when a is undefined or 0-by-0
                if all(sza==[0 0])
                    sza(2) = 1;
                end
                L(1).subs{2} = 1:sza(2);
            end
            b = PVsubsasgn_2idx(a,L,temp);
        end
        
    case '{}'
        error('{}- like subsassign is not supported for polynomial objects.');
end
