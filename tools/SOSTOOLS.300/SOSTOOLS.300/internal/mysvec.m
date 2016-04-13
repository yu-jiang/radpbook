function x = mysvec(X)
% This function is my implementation of the svec function of SDPT3

% SP

s = size(X,1);
ii = ones(s);
idx1 = find(triu(ii));
idx2 = find(triu(ii,1));
X(idx2) = X(idx2)*sqrt(2);
x = X(idx1);