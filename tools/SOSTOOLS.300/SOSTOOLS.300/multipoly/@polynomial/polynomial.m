% function P = polynomial(Coefficient,Degmat,Varname,Matdim);
%
% DESCRIPTION
%   Creates a polynomial or a matrix of polynomials.
%
% INPUTS
%   Coefficient: coefficients of each monomial.
%   Degmat: degrees of each monomial
%   Varname: names of variables
%   Matdim: dimensions of the polynomial matrix
%
% OUTPUT
%   P: polynomial object
%
% SYNTAX
%   P=polynomial
%     Creates an empty polynomial object.
%   P=polynomial(Coefficient)
%     If Coefficient is a real matrix of dimension NxM, then P is
%     an NxM constant polynomial.
%   P=polynomial(Varname)
%     If Var is an NxM cell array of strings, then P is an NxM polynomial
%     whose entries are the variables specified in Var.
%   P=polynomial(Coefficient)
%     If Coefficient is a polynomial object, then P=Coefficient.
%   P=polynomial(Coefficient,Degmat,Varname,Matdim)
%     If P is an NxM polynomial that is the sum of T terms in V
%     variables, the inputs should be specified as:
%       Coefficient is a Tx(N*M) sparse matrix.
%          The coefficients of the (i,j) entry of the polynomial
%          matrix are a Tx1 vector stored in the i+*N*(j-1) column of
%          Coefficient.
%       Degmat is a TxV sparse matrix of natural numbers.  Row t
%          gives the degrees of each variable for the t^th term.
%       Varname is a Vx1 cell array with entry v giving the name of
%          variable v. For a constant polynomial, varname is an
%          empty 1x1 cell.
%       Matdim is a 1x2 vector of the matrix dimensions, [N M].

% 6/7/2002:   PJS  Initial Coding
% 6/8/2002:   PJS  Allow Matrices of Polynomials
% 10/30/2002: PJS  Input is structure not cell arrays
% 11/5/2002:  PJS  Coefficients stored as sparse 3D array.
% 4/20/2009:  PJS  Changed Coefficient matrix from (T*N*M)x1 to Tx(N*M)
% 5/20/2009:  PJS  Added chkval (=0 to skip validity check, else = 1)
% 10/21/2010: PJS  Converted to new Matlab format for object class defs


classdef polynomial
    % Properties
    properties
        coefficient = sparse([]);
        degmat = sparse([]);
        varname = {};
    end
    
    % XXX PJS 10/21/2010:
    % I need to create my own subsref.  My subsref is a method and hence
    % it has access to private properties. Since '.'-refs go through
    % subsref, I need to handle priviledge checking in my subsref. For
    % now I'm sticking with AP's get/set functionality.
    properties (SetAccess = private)
        % XXX Change PVget/PVset so that matdim is not get/setable?
        matdim = [0 0];
    end
    properties (SetAccess = private, Dependent = true)
        nterms;
        nvars;
        maxdeg;
        mindeg;
    end
    
    % Methods:
    % All methods except the constructor are defined in separate files
    % Note: The private folder in @polynomial is not needed. I could
    % specify the methods as private in this classdef file.
    methods
        % Constructor
        function p = polynomial(coefficient,degmat,varname,matdim,chkval)
            %--------------------------------------%
            %    Fill in empty arguments
            %--------------------------------------%
            if nargin>=1 && isa(coefficient,'logical');
                % Lift boolean to a double
                coefficient = double(coefficient);
            end
            
            if nargin == 0
                % Return p with default values
                return
            elseif nargin==1
                if isa(coefficient,'polynomial')
                    % If input is already a polynomial, pass it back
                    p = coefficient;
                    return;
                elseif isa(coefficient,'double') && isreal(coefficient) && ...
                        ndims(coefficient)==2
                    
                    % If input is a real double, then convert to a poly
                    szc = size(coefficient);
                    p.coefficient = reshape(coefficient,[1 szc(1)*szc(2)]);
                    p.degmat = sparse(1,0);
                    p.varname  = {};
                    p.matdim = szc;
                    return
                elseif iscellstr(coefficient) && ndims(coefficient)==2
                    V = coefficient;
                    szV = size(V);
                    numV = szV(1)*szV(2);
                    
                    p.coefficient = speye(numV);
                    p.degmat = speye(numV);
                    p.varname = V(:);
                    p.matdim = szV;
                    return;
                else
                    error(['For single input, argument must be a polynomial, '...
                        'real double, or cell array of strings']);
                end
            elseif nargin<4
                errstr1 = 'Invalid syntax for the "polynomial" command.';
                errstr2 = ' Type "help polynomial" for more information.';
                error([errstr1 errstr2]);
            elseif nargin==4
                chkval=1;
            end
            
            %--------------------------------------%
            %       Error Checking
            %--------------------------------------%
            p.coefficient = coefficient;
            p.degmat = degmat;
            p.varname = varname(:);
            p.matdim = matdim(:)';
            
            if chkval
                [flag,errormsg] = PVisvalid(p);
            else
                % Skip validity check -- This is used by internal
                % functions, e.g. combine, for which validity check is not
                % needed.  This speeds up the call to polynomial.
                flag = 1;
            end
            
            if flag==1
                % p is a valid polynomial object
                p.coefficient = sparse(p.coefficient);
                p.degmat = sparse(p.degmat);
            else
                % p is not a valid polynomial object.
                error(errormsg);
            end
        end % polynomial constructor
        
    end % methods
end % classdef






