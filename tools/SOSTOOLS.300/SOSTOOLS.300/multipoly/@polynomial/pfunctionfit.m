function [pfit,cfit,info] = pfunctionfit(p,vars,Xdata,fnc,W)
% function [pfit,cfit,fiterr] = pfunctionfit(p,x,Xdata,fnc,W)
%
% DESCRIPTION:
%   This function finds the coefficients of a multivariate polynomial that
%   best fits a function fnc in a least-squares cost.  The function is fit
%   with a linear combination of polynomial basis functions:
%      p(x,c) = c1*f1(x)+c2*f2(x) + ... + ck*fk(x)
%   where f1, f2, ..., fk are the polynomial basis functions. pfunctionfit
%   samples the function fnc and computes the coefficients c1, c2, ..., ck
%   that minimize the fitting error on these samples in a weighted squares
%   squares cost.  See pdatafit for more detail.
%
% INPUTS:
%   p:     1-by-1 polynomial.
%   x:     Nx-by-1 vector of pvars that specifies the independent variables
%          in p.  All other variables in p are considered to be coefficients
%   Xdata  (Optional): Nx-by-Npts  matrix of input data values at which to
%          evaluate fnc for fitting. Alternatively, Xdata can be a 
%          structure with fields specying how to construct the data:
%          fields to construct the Nvars-by-Npts input data set.
%           - range:= Nx-by-2 matrix containing the data range [min max] of
%             the variables defined in vars. Default is [-1 1]
%           - sample:=defines the sampling technique. Choices are: 'grid', 
%             'uniform' ,'lhs'. Default is 'grid'. 'grid' generates 
%             linearly spaced data along each direction. 'uniform' draws
%             random samples from the range using a uniform distribution.
%             'lhs' uses the Latin Hypercube sampling technique. 'lhs'
%              requires the Statistics Toolbox.
%           - Npts: If sample='grid' then Npts is a Nvar-by-1 vector 
%             defining the number of points to be sampled along each
%             direction. The total # of points is prod(Npts).  For 'lhs'
%             or 'uniform', Npts is a 1-by-1 defining the total
%             number of sampled points.
%  fnc :   Function to fit with inputs x and 1-by-1 output. fnc can be
%          a function handle, string expression, or polynomial.
%  W (Optional): 1-by-Npts weighting vector .  Alternatively W can be
%          a function handle, string expression, or polynomial.
%          [Default: W=ones(1,Npts)]
%
% OUTPUTS:
%   pfit: Least-squares polynomial fit
%   cfit: Nc-by-2 cell array of the optimal coefficients.  The first
%         column contains the coefficients (as chars) and the second
%         column contains the optimal values.  The subs command can be
%         be used to replace the coefficients in any polynomial with
%         their optimal values, e.g. pfit = subs(p,cfit).
%   info: Data structure containing the matrices in the least squares
%         problem.  info has the fields A, b, cfit, W, e as described
%         in pdatatfit help. It also contains Xdata and Ydata. Xdata are
%         the input data samples and Ydata gives the values of fnc
%         evaluated at Xdata. Sample information is stored in the fields
%         sample, range, and Npts.
%
% SYNTAX
%   [pfit,cfit,info] = pfunctionfit(p,vars,Xdata,fnc)
%   [pfit,cfit,info] = pfunctionfit(p,vars,fnc)
%   [pfit,cfit,info] = pfunctionfit(p,vars,Xdata,fnc,W)
%
% EXAMPLE
%  fnc = @(x) sin(x);
%  pvar c0 c1 c2 c3 x;
%  vars = x;
%  p = c0 + c1*x + c2*x^2 + c3*x^3;
%  Xdata.sample ='uniform';
%  Xdata.Npts = 20;
%  [pfit,cfit,info] = pfunctionfit(p,vars,Xdata,fnc);
%  ezplot(fnc,[info.range(1) info.range(2)]); hold on;
%  xx = linspace(info.range(1),info.range(2),10);
%  plot(xx,double(subs(pfit,vars,xx)),'r--'); hold off;
%  legend('Original Function','Polynomial Fit'); xlabel('x');
%
% See also pdatafit, lhsdesign


% Abhijit  12/07/09   Initial Coding
% PJS      11/23/10   Updates to help and minor mods to code


% ------------ Error Checking

% ------ Input Argument Handle
if nargin == 1 || nargin == 2
    error('Invalid syntax for pfunctionfit. Type "help pfunctionfit" for more information.');
elseif nargin == 3
    fnc = Xdata; Xdata = []; W = [];
elseif nargin == 4
    W = [];
end

% ------- Checking if the function is a handle or polynomial object
if ~( isa(fnc,'function_handle') || isa(p,'polynomial') || ischar(fnc) )
    error('Function must be a polynomial, function handle, or string')
end


% ---------- Error Check if fitting is to be done with undefined variables
if isa(fnc,'polynomial') && ~isempty(setdiff(char(vars),fnc.varname))
    error('Variables in fnc are inconsistent with those listed in var');
end


% --------- Handle the Function Format
% Handle string function: Posible two cases  (i) built-in function
% ('sin'...)  (ii) string like '3*x + 5'
if ischar(fnc)
    if exist(fnc,'builtin')== 5  % built-in
        fnc = str2func(fnc);
    else
        fnc = inline(fnc); % string like (ii)
    end
end

% ---- handle W if it is a function
if ~isempty(W) && ischar(W)
    if exist(W,'builtin')== 5  % built-in
        W = str2func(W);
    else
        W = inline(W); % string like (ii)
    end
end

% ---- Size of variable
Nvar = length(vars);


% ------------ Handle Xdata : Could be either data set or structure to
% define how to construct dataset.


% --- Handle if Xdata is data provided by user
if ~isstruct(Xdata) && ~isempty(Xdata)
    if size(Xdata,1) ~= Nvar
        error('Xdata should be Nvar-by-Npts')
    end
    % ---- W can only be provided if Xdata is provided
    % Can't W be a function handle?
    %if isempty(Xdata) && ~isempty(W)
    %    error('Weighting can only be provided if Xdata is provided')
    %end
    
elseif isstruct(Xdata) || isempty(Xdata)
    % ---- check if its a structure or not and fill out information
    
    % ---- Default Sampling Field
    if ~isfield(Xdata,'sample') || isempty(Xdata.sample)
        Xdata.sample = 'grid';
    end
    
    % ---- Default Data Range
    if ~isfield(Xdata,'range') || isempty(Xdata.range)
        Xdata.range = repmat([-1 1],Nvar,1);
    end
    
    % -------- is enough information given about range
    if Nvar ~= size(Xdata.range,1)
        error('Dimension mismatch between vars and range')
    end
    
    % --------- range is not well defined
    if any( Xdata.range(:,1) > Xdata.range(:,2) )
        error('Data range is not well defined')
    end
    
    % ---- Default Npts
    if ~isfield(Xdata,'Npts') || isempty(Xdata.Npts)
        
        if strcmp(Xdata.sample,'uniform') || strcmp(Xdata.sample,'lhs')
            % --- Uniform and LHS can have 1-by-1 total no of datapoints
            Xdata.Npts = 10^Nvar;
            
        elseif strcmp(Xdata.sample,'grid')
            % ---- grid sampling specifies an array of Npts for each dimensionn
            Xdata.Npts = repmat(10,Nvar,1);  % 10 datapoints along each direction
        end
        
    elseif isfield(Xdata,'Npts') || ~isempty(Xdata.Npts)
        if strcmp(Xdata.sample,'uniform') || strcmp(Xdata.sample,'lhs')
            % --- Uniform and LHS can have 1-by-1 total no of datapoints
            if any( size(Xdata.Npts) ~= [1 1] )
                error(' For uniform and LHS sampling only provide total no of points')
            end
        elseif strcmp(Xdata.sample,'grid')
            % ---- grid sampling specifies an array of Npts for each dimensionn
            if any(size(Xdata.Npts)~=1)  && any(size(Xdata.Npts)~=Nvar)
                error('Incorrect dimension of NPTS for GRID Sampling')
            end
        end
    end
    
end


%---------- Handle the Weighting Function
if ~isempty(W)
    %  Check if Weighting function is to be done with undefined variables
    if  isa(W,'polynomial') && ~isempty(setdiff(W.varname, char(vars)))
        error('Weighting Function contains undefined variables')
    end
end

%------------------------------------------- Generate Data Set
if isstruct(Xdata) && strcmp(Xdata.sample,'grid')
    % -------------- Gridded data set
    range = Xdata.range;
    Npts = Xdata.Npts
    if Nvar == 1
        Xdata.X = linspace(range(1,1),range(1,2),Npts);
    else
        X = cell(1,Nvar);
        for i1 = 1:Nvar
            X{i1} = linspace(range(i1,1),range(i1,2),Npts(i1));
        end
        
        % --- Generate Gridded Data
        Xgrid = cell(1,Nvar);
        [Xgrid{:}]=ndgrid(X{:});
        Xdata.X = reshape(cat(length(Xgrid)+1,Xgrid{:}),[],length(Xgrid)) ;
        Xdata.X = Xdata.X';
        
        % -- total npts for grid
        Npts = prod(Xdata.Npts);
        Xdata.Npts = Npts;
    end
    
elseif isstruct(Xdata) && strcmp(Xdata.sample,'uniform')
    % --------------- Uniform data set
    range = Xdata.range;
    Npts = Xdata.Npts;
    if Nvar == 1
        Xdata.X = range(1,1) + (range(1,2) - range(1,1))*rand(1,Npts);
    else
        % creating data within the range... x = a + (b-a)*samplepoint
        AA = repmat(range(:,1),1,Npts);  % the min range of all variables
        BB = repmat( (range(:,2) - range(:,1)),1,Npts);
        Xdata.X = AA + BB.*rand(Nvar,Npts);
    end
        
elseif ~isstruct(Xdata) && strcmp(class(Xdata),'double')
    % -------------- User Supplied data set
    
    if any(size(Xdata)==[1 1])
        % Make sure Xdata lines up for single variable case
        % (If there are multiple variables then user must make sure
        %  dimensions of Xdata are as described in the help)
        Xdata.X= Xdata(:)';
        Xdata.Npts = length(Xdata.X);
    else
        Xdata.Npts = size(Xdata,2);        
        if size(Xdata,1) ~= Nvar
            error('Xdata should be dimension of Nvar-by-Npts')
        else
            Xdata.X = Xdata;
        end
    end
    Npts = Xdata.Npts;
    
elseif isstruct(Xdata) &&  strcmp(Xdata.sample,'lhs')
    % -------------- Latin HyperCube Sampling
    range = Xdata.range;
    Npts = Xdata.Npts;
    
    % Check if lhsdesign is installed in matlab or not...
    if exist('lhsdesign','file')~= 0
        if Nvar == 1
            Xdata.X = range(1,1) + (range(1,2) - range(1,1))*lhsdesign(1,Npts);
        else
            % creating data within the range... x = a + (b-a)*samplepoint
            AA = repmat(range(:,1),1,Npts);  % the min range of all variables
            BB = repmat( (range(:,2) - range(:,1)),1,Npts);
            Xdata.X = AA + BB.*lhsdesign(Nvar, Npts);
        end        
    else
        error('STATISTICS TOOLBOX is required to use LHS sampling')
    end
end

% --- Generate Weighting Array
if ~isempty(W)
    if isa(W,'polynomial')
        W = double(subs(W,vars,Xdata.X));
    elseif isa(W,'function_handle') || strcmp(class(W),'inline')
        W =  LOCALWeval(W,Xdata.X,Xdata.Npts);
    elseif isa(W,'double')
        W = W(:);
        if length(W)~=Npts
            error('Weighting Dimension mismatch with provided data')
        end
    end
else
    W = ones(Xdata.Npts,1);
end

%------------------------ Generate Ydata
if isa(fnc,'polynomial')
    Ydata = double(subs(fnc,vars,Xdata.X));
elseif isa(fnc,'function_handle') || isa(fnc,'inline')    
    Ydata = zeros(1,Xdata.Npts);    
    if nargin(fnc) == 1
        for i1 = 1:Xdata.Npts
            Ydata(i1) = fnc(Xdata.X(:,i1));
        end        
    else
        % Convert all data to cell array
        % Xdatacell = mat2cell(Xdata,repmat(1,1,Nvar),repmat(1,1,Npts));
        Xdatacell = num2cell(Xdata.X);
        for i1 = 1:Xdata.Npts
            Ydata(i1) = fnc(Xdatacell{:,i1});
        end        
    end
end

%------------------------- Use Pdatafit
[pfit,cfit,info] = pdatafit(p,vars,Xdata.X,Ydata,W);
info.Ydata = Ydata;
info.Xdata = Xdata.X;  
info.range =[]; 
if isfield(Xdata,'range') 
    info.range = Xdata.range;
end
info.Npts =[]; 
if isfield(Xdata,'Npts')
    info.Npts = Xdata.Npts; 
end
info.sample =[]; 
if isfield(Xdata,'sample') 
    info.sample = Xdata.sample; 
end

% ------ Local Function to Evaluate W if it is a function handle
function [Wdata] =  LOCALWeval(W,XX,XNpts)

Wdata = zeros(1,XNpts);
if nargin(W) == 1    
    for i1 = 1:XNpts
        Wdata(i1) = W(XX(:,i1));
    end    
else
    % Convert all data to cell array
    Xdatacell = num2cell(XX);
    for i1 = 1:XNpts
        Wdata(i1) = W(Xdatacell{:,i1});
    end    
end