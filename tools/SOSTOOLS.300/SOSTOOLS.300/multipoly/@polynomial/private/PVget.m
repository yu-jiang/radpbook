function out  = PVget(a,property)
% function Out = PVget(A,Property)
%
% DESCRIPTION
%   (INTERNAL FUNCTION)
%   Private get for polynomial objects.
%
% INPUTS
%   A: a polynomial
%   Property: polynomial property
%
% OUTPUTS
%   Out: Value of the property
%
% SYNTAX
%   Out = PVget(A,Property);

% 6/7/2002: PJS  Initial Coding
% 10/31/2002 PJS Added properties for scalar polynomials
% 1/31/2003: SP  Added 'matdim' as a getable/setable property
% 3/30/2003: PJS Removed 'matdim' as a getable/setable property
%                (use size and reshape instead)
% 2/27/2008: PJS  Properly order varname output for vector of pvars
% 11/10/2008: PJS  Undo the 2/27/2008 update because of some bugs

switch property
    
    case 'PropNames'
        % These are getable and setable PUBLIC properties
        out.GPropNames = {'coefficient' ; 'degmat';  'varname'; ...
            'nterms'; 'nvars'; 'matdim'; 'maxdeg'; 'mindeg'};
        
        out.SPropNames = {'coefficient' ; 'degmat'; 'varname'};
        
    case 'coefficient'
        out = a.coefficient;
        
    case 'degmat'
        out = a.degmat;
        
    case 'varname'
        out = a.varname;
        
    case 'nterms'
        out = size(a.degmat,1);
        
    case 'nvars'
        out = size(a.degmat,2);
        
    case 'maxdeg'
        if all( size(a.degmat)==[0 0]) %isempty(a.degmat)
            out = [];
        else
            out = full(max(sum(a.degmat,2)));
        end
        
    case 'mindeg'
        if all( size(a.degmat)==[0 0]) %isempty(a.degmat)
            out = [];
        else
            out = full(min(sum(a.degmat,2)));
        end
        
    case 'matdim'
        out = a.matdim;
        
    otherwise
        error(['Property ' property ' is Invalid']);
        
end



