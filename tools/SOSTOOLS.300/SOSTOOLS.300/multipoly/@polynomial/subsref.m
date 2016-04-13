function varargout = subsref(a,L)
% function B=subsref(A,L);
%
% DESCRIPTION
%   Subsreference for polynomial objects.
%
% INPUTS
%   A: polynomial
%   L: a structure array with the fields:
%    type -- string containing '()'  or '.' specifying the
%              subscript type.
%    subs -- Cell array or string containing the actual subscripts.
%
% OUTPUTS
%   B: object after referencing
%
% SYNTAX
%   B=subsref(A,L)

% 6/7/2002: PJS  Initial Coding

switch L(1).type
    case '.'
        b = get(a,L(1).subs);
    case '()'
        sza = a.matdim;
        nra = sza(1);
        nca = sza(2);
        acoef = a.coefficient;
        b = a;
        
        if length(L(1).subs)==1
            % Check Indices
            subsidx = L(1).subs{1};
                
            if strcmp(subsidx,':')
                subsidx = (1:nra*nca)';
            end
            if islogical(subsidx)
                subsidx = find(subsidx);
            end
            
            if ~isempty(subsidx) && ...
                    (max(subsidx)>nra*nca || min(subsidx)<1)
                error('Index exceeds matrix dimensions.');
            end
            
            % Do subsref
            b.coefficient = acoef(:,subsidx);
            if isempty(subsidx)
                b = polynomial(zeros(size(subsidx)));                
            elseif nca==1
                % Column vector should remain column vector
                b.matdim = [length(subsidx) 1];
            else
                % All other single index references convert to row vector
                b.matdim = size(subsidx);
            end
    
        elseif length(L(1).subs)==2
            ridx = L(1).subs{1};
            if strcmp(ridx,':')
                ridx = 1:nra;
            end
            cidx = L(1).subs{2};
            if strcmp(cidx,':');
                cidx = 1:nca;
            end
            
            % Check Indices
            if ~isempty(ridx) && (max(ridx)>nra || min(ridx)<1)
                error('Index exceeds matrix dimensions.');
            elseif ~isempty(cidx) && (max(cidx)>nca || min(cidx)<1)
                error('Index exceeds matrix dimensions.');
            end
            
            % Do subsref
            bcoef = acoef;
            idx = [];
            for i1 = 1:length(cidx)
                idx = [idx ridx+(cidx(i1)-1)*nra];
            end
            b.coefficient = bcoef(:,idx);
            b.matdim = [length(ridx) length(cidx)];
            
        else
            error('Invalid subsref for polynomials');
        end
        
%         % Remove any terms with zero coeffs and then combine down
%         % (This is faster than directly calling combine without
%         %  pruning out terms with zero coefs)
%         bcoef = b.coefficient;
%         [ridx,cidx] = find( bcoef );
%         ridx = unique(ridx);
%         b.coefficient = bcoef(ridx,:);
%         b.degmat = b.degmat(ridx,:);
        
        if isempty(b.coefficient)
            b = polynomial(zeros(b.matdim));
        else
            b = combine(b);
        end
        
    case '{}'
        error('{}- like subsreference is not supported for polynomial objects.');
end

if length(L)>1
    b = subsref(b,L(2:end));
end
varargout{1}=b;


