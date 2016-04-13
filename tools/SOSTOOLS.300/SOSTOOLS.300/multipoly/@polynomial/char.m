function sr = char(a)
% function SR=char(A)
%
% DESCRIPTION
%   Converts a polynomial to its string representation.
%
% INPUTS
%   A: r-by-c polynomial
%
% OUTPUTS
%   SR: r-by-c cell array with each entry giving the string representation
%       of the corresponding entry of the polynomial.
%
% SYNTAX
%   SR=char(A)

% 6/9/2002 PJS Initial Coding

% Create String Representation
sza = size(a);
if min(sza)==0
    sr = {['Empty polynomial: ' int2str(sza(1)) '-by-' int2str(sza(2))]};
    return;
end

% Fast method to convert arrays of pvars
if ispvar(a)
    v = a.varname;
    [ridx,cidx]=find(a.degmat);
    v(ridx) = v(cidx);
    
    acoef = a.coefficient;
    [ridx,cidx]=find(acoef);
    v(cidx) = v(ridx);
    
    sr = reshape(v,sza);
    return;
end

% Get info from polynomial structure
acoef = a.coefficient;
adeg = a.degmat;
avar = a.varname;
nt = size(adeg,1);
nv = size(adeg,2);

sr = cell(sza);
for i1 = 1:sza(1);
    for i2 = 1:sza(2);
        
        % Make a char from the (i1,i2) entry
        s = [];
        firstterm = 1;
        idx = nt*(i1-1) + nt*sza(1)*(i2-1) + (1:nt);
        for i3 = 1:nt
            
            if acoef(idx(i3))~=0
                % Sign of term
                if firstterm~=1 && acoef(idx(i3))>=0
                    s = [s ' + '];
                elseif firstterm~=1 && acoef(idx(i3))<0
                    s = [s ' - '];
                elseif firstterm==1 && acoef(idx(i3))<0
                    s = '-';
                end
                firstterm=0;
                
                if nv==0 || sum(adeg(i3,:))==0
                    % constant term
                    s = [s num2str(abs(acoef(idx(i3))))];
                else
                    % not constant term
                    if abs(acoef(idx(i3))) ~=1
                        s = [s num2str(abs( acoef(idx(i3)) )) '*'];
                    end
                    for i4 = 1:nv
                        if adeg(i3,i4)>1
                            % XXX R2010b -- int2str does not work for
                            % sparse inputs.
                            s = [s avar{i4} '^' int2str( full(adeg(i3,i4)) ) '*'];
                        elseif adeg(i3,i4)==1
                            s = [s avar{i4} '*'];
                        end
                        if i4==nv
                            s = s(1:end-1);
                        end
                    end
                end
            end
            
        end
        
        if firstterm==0
            sr{i1,i2} = s;
        else
            sr{i1,i2} = '0';
        end
        
    end
end




