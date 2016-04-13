function [A,b,K,z,dv2x,Nfv,feas,zrem] = sospsimplify(A,b,K,z,dv2x,Nsosvarc)
% function [A,b,K,z,dv2x,Nfv,feas,zrem] = sospsimplify(A,b,K,z,dv2x,Nsosvarc)
%
% DESCRIPTION 
%   This function performs a simplification procedure on the SOS problem.
%   First, it tries to detect the sign of optimization variables based on
%   simple constraints. Second, it searches for monomials that can be 
%   removed from each SOS constraint. This search is based on diagonal 
%   entries of the Gram matrix that are forced to be zero and it is 
%   equivalent to the Newton Polytope method. These two steps are repeated
%   until no new sign information can be detected.  The removed monomials
%   are stored in zrem.
%
%   In the code, the information about the sign of the optimization 
%   variables is stored in xsign where:
%      xsign(i)=NaN    if x(i) has unknown sign
%      xsign(i)=+1     if x(i)>=0
%      xsign(i)=-1     if x(i)<=0
%      xsign(i)=0      if x==0
%

% 9/21/08  PJS Initial Coding
% 12/14/10 PJS Bug related to re-indexing of dv2x when removing free vars  


%--------------------------------------------------------------------
% Grab problem dimensions
%--------------------------------------------------------------------
Nfv = K.f;
Nlp = K.l;
%Nsv = sum( K.s(1:Nsosvarc).^2 );
Nx = size(A,2);

%--------------------------------------------------------------------
% Find non-negative optimization variables
%--------------------------------------------------------------------
xsign = NaN([Nx,1]);

% Process LP constraints: A*d+y=b, y>=0
if Nlp>0
    % LP slack vars are >= 0
    xsign(Nfv+1:Nfv+Nlp) = +1;
end

% Process SOS Inequality Constraints: Ad*d + Aq*Q(:) = b and Q>=0
ptr = Nfv+Nlp;
for i1=1:length(K.s)
    % Diag entries of Q are >= 0    
    lz = K.s(i1);
    diagidx = (0:lz-1)*lz+(1:lz);
    xsign(ptr+diagidx) = +1;
    ptr = ptr+lz^2;
end

%--------------------------------------------------------------------
% Find dec vars and monomials that can be removed
%--------------------------------------------------------------------
go = 1;
zrem = cell(length(K.s),1);
xsignprev = xsign;
while (go == 1)
    % Use simple constraints to determine sign of optimization vars. 
    xsign = LOCALxsignupdate(xsign,A,b);
    A( : , xsign==0 ) = 0;
    
    % Monomial reduction procedure 
    % (This is equivalent to the Newton Polytope method)
    ptr = Nfv+Nlp;
    for i1=1:length(K.s)
        % Find diag entries of Q that are forced to be zero
        lz = K.s(i1);
        blkidx = ptr+(1:lz^2);
        Qsign = diag( reshape(xsign(blkidx),[lz lz]) );
        loc = find(Qsign==0);        
        if ~isempty(loc)
            % Corresponding rows/cols of zero diag entries are also zero
            tmp = sparse(lz,lz);
            tmp(loc,:) = 1;
            tmp(:,loc) = 1;
            rmidx = find(tmp)+ptr;
            xsign(rmidx) = 0;    
            
            % Remove vars/monoms associated with zero Gram matrix entries
            A(:,rmidx) = [];
            xsign(rmidx) = [];
            zrem{i1} = [zrem{i1}; z{i1}(loc)];
            z{i1}(loc) = [];
            %K.s(i1) = size( z{i1}, 1); %length( z{i1} );
            if isempty(z{i1})
                K.s(i1) = 0;
            else
                K.s(i1) = length( z{i1} );
            end
            
            % Update the mapping of dec vars into the optimization vars
            if i1<=Nsosvarc
                % Number of optim vars to remove
                Nremove = length(rmidx);
                
                % Optim vars currently in this block (before removal)
                blkidx = ptr+(1:lz^2);

                % Mark removed dec vars
                dv2x( ismember(dv2x,rmidx) ) = 0;
                          
                % Relabel the remaining dec vars in this block
                idx = find( tril(ones(K.s(i1))) );                            
                dv2x( ismember(dv2x,blkidx) ) = ptr+idx;
                                
                % Relabel remaining dec vars in subsequent blocks
                idx = find( dv2x>blkidx(end) );
                dv2x(idx) = dv2x(idx) - Nremove;                                               
            end        
        end
        
        % Update pointer
        ptr = ptr+K.s(i1)^2;
    end

    % Continue if xsign has been updated
    go = ~isequalwithequalnans(xsign,xsignprev);
    xsignprev = xsign;    
end

%--------------------------------------------------------------------
% Clean up 
%--------------------------------------------------------------------

% Mark removed free decision vars
rmidx = find( xsign(1:K.f)==0 );
dv2x( ismember(dv2x,rmidx) ) = 0;
idx = find(dv2x<=K.f & dv2x>0);
A(:,rmidx) = [];
xsign(rmidx) = [];
Nremf = length(rmidx);
K.f = K.f - Nremf;  
%K.f = K.f - length(rmidx);  
Nfv = K.f;
dv2x(idx) = 1:Nfv;

idx2 = find(dv2x>K.f & dv2x>0);
dv2x(idx2) = dv2x(idx2)-Nremf;

% Remove any constraints of the form 0=0 
% XXX -- We need a tolerance here.  How should we choose tol?
tol = 1e-9;
ridx = find( sum(A~=0,2)==0 & abs(b)<max(tol,tol*max(abs(b))) );
A(ridx,:) = [];
b(ridx) = [];

% Check for infeasible problems of the form 0 = bi where bi is not equal 
% to zero (Our simplify code should flag infeasible problems because
% Sedumi can error out on problems that are trivially infeasible)
if isempty(A)
    feas = 1;
    return    
else
    ridx = find( sum(A~=0,2)==0 & abs(b)>tol*max(abs(b)) );
end
feas = 1;
if ~isempty(ridx)
    feas = 0;
end

% Add equality constraints for remaining dec vars known to be zero
idx = find( xsign==0 );
lidx = length(idx);
A(end+1:end+lidx,idx) = speye(lidx);
b = [b; sparse(lidx,1)];



%--------------------------------------------------------------------
% Local function to update sign of optimization var 
%--------------------------------------------------------------------
function xsign = LOCALxsignupdate(xsignOld,A,b)

% Initialize output
xsign = xsignOld;

% Process constraints of the form:  aij*xj = bi
ridx = find(  sum(A~=0,2)==1  );
if ~isempty(ridx)
    [cidx,tmp]=find( A(ridx,:)' );
    idx = sub2ind(size(A),ridx,cidx);

    signA = sign( A(idx) );
    signb = sign( b(ridx) );    
    xsignUpdate = signA.*signb;
    
    % XXX PJS 12/07/09: If cidx = [2;2] and xsignUpdate = [1; NaN] then 
    % the next line will replace xsign(2) with NaN because the last index 
    % in a subsasgn wins. This caused problems on a GSOSOPT problem.
    %
    %xsign(cidx) = LOCALupdate(xsign(cidx),xsignUpdate);
    
    % The correct code (also below) is below.  I'll try to vectorize
    % if speed becomes an issue.
    for i1 =1:length(cidx)
        xsign(cidx(i1)) = LOCALupdate(xsign(cidx(i1)),xsignUpdate(i1));
    end
        
end

% Process constraints of the form:  aij*xj + aik*xk = bi
ridx = find(  sum(A~=0,2)==2  );
if ~isempty(ridx)
    [cidx,tmp]=find( A(ridx,:)' );
    cidx = reshape(cidx,[2 length(ridx)])';
    cidx1 = cidx(:,1);
    idx1 = sub2ind(size(A),ridx,cidx1);
    cidx2 = cidx(:,2);
    idx2 = sub2ind(size(A),ridx,cidx2);

    c1 = b(ridx)./A(idx1);
    c2 = -A(idx2)./A(idx1);
    xsignUpdate = NaN([length(ridx) 1]);
    xsignUpdate( c1<=0 & (c2.*xsign(cidx2)<=0) ) = -1;
    xsignUpdate( c1>=0 & (c2.*xsign(cidx2)>=0) ) = +1;
    %xsign(cidx1) = LOCALupdate(xsign(cidx1),xsignUpdate);
    for i1 =1:length(cidx1)
        xsign(cidx1(i1)) = LOCALupdate(xsign(cidx1(i1)),xsignUpdate(i1));
    end
    
    c1 = b(ridx)./A(idx2);
    c2 = -A(idx1)./A(idx2);
    xsignUpdate = NaN([length(ridx) 1]);
    xsignUpdate( c1<=0 & (c2.*xsign(cidx1)<=0) ) = -1;
    xsignUpdate( c1>=0 & (c2.*xsign(cidx1)>=0) ) = +1;
    %xsign(cidx2) = LOCALupdate(xsign(cidx2),xsignUpdate);
    for i1 =1:length(cidx2)
        xsign(cidx2(i1)) = LOCALupdate(xsign(cidx2(i1)),xsignUpdate(i1));
    end        
end

% Process constraints of the form:  aij*xj + aik*xk + ail*xl= 0
% where aij*xj, aik*xk, ail*xl all have the same sign.  
% This implies that each of the three vars  = 0
ridx = find(  sum(A~=0,2)==3 & b==0 );
if ~isempty(ridx)
    [cidx,tmp]=find( A(ridx,:)' );
    cidx = reshape(cidx,[3 length(ridx)])';
    cidx1 = cidx(:,1);
    idx1 = sub2ind(size(A),ridx,cidx1);
    cidx2 = cidx(:,2);
    idx2 = sub2ind(size(A),ridx,cidx2);
    cidx3 = cidx(:,3);
    idx3 = sub2ind(size(A),ridx,cidx3);

    % All terms are non-neg
    rsign = (A(idx1).*xsign(cidx1)>=0) & (A(idx2).*xsign(cidx2)>=0) ...
                & (A(idx3).*xsign(cidx3)>=0);
    idx = find(rsign==1);    
    for i1=idx
        xsign(cidx1(i1)) = 0;
        xsign(cidx2(i1)) = 0;
        xsign(cidx3(i1)) = 0;        
    end
    
    % All terms are non-pos
    rsign = (A(idx1).*xsign(cidx1)<=0) & (A(idx2).*xsign(cidx2)<=0) ...
                & (A(idx3).*xsign(cidx3)<=0);
    idx = find(rsign==1);    
    for i1=idx
        xsign(cidx1(i1)) = 0;
        xsign(cidx2(i1)) = 0;
        xsign(cidx3(i1)) = 0;        
    end
end

%--------------------------------------------------------------------
% Local function to update sign of optimization var
%--------------------------------------------------------------------
function xsignNew = LOCALupdate(xsign,xsignUpdate)

% Find constraints that force ai=0, ai<0 and ai>0
zidx = find( (xsignUpdate==0) | (xsign==-1 & xsignUpdate>0) | ...
    (xsign==+1 & xsignUpdate<0) );
nidx = find( isnan(xsign) & xsignUpdate<0 );
pidx = find( isnan(xsign) & xsignUpdate>0 );

% Update xsign
xsignNew = xsign;
xsignNew( zidx ) = 0;
xsignNew( nidx ) = -1;
xsignNew( pidx ) = +1;

