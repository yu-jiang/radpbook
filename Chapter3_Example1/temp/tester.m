% Test if the A B Matrices are valid

% Coefficients
mdlCoeff


A = [ 0 1 0 0; [-ks -bs ks bs]/mb ; ...
      0 0 0 1; [ks bs -ks-kt -bs]/mw];
B = [0; 10000/mb ; 0;  -10000/mw];

susp_sys_linear(rand(4,1),rand, A, B);

P_init = lyap(A',eye(4));
K1 = 1/r*B'*P_init;

[P,~,K] = care(A,B,eye(4),1)



% Optimal values 
%
% P =
% 
%     1.8639    0.0703   -1.5871    0.0123
%     0.0703    0.0638   -0.9956    0.0069
%    -1.5871   -0.9956   40.3489   -0.0763
%     0.0123    0.0069   -0.0763    0.0063
% 
% 
% K =
% 
%     0.2868    0.9727  -20.4658   -0.8259