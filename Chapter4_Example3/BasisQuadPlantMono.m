function y = BasisQuadPlantMono(x1,x2)
% Basis for the generalized HJ
y = [x1*x1;
    x1*x2;
    x2*x2;
    x1*x1*x1;
    x1*x1*x2;
    x1*x2*x2;
    x2*x2*x2;
    x1*x1*x1*x1;
    x1*x1*x1*x2;
    x1*x1*x2*x2;
    x1*x2*x2*x2;
    x2*x2*x2*x2;
    x1*x1*x1*x1*x1;
    x1*x1*x1*x1*x2;
    x1*x1*x1*x2*x2;
    x1*x1*x2*x2*x2;
    x1*x2*x2*x2*x2;
    x2*x2*x2*x2*x2;
    x1*x1*x1*x1*x1*x1;
    x1*x1*x1*x1*x1*x2;
    x1*x1*x1*x1*x2*x2;
    x1*x1*x1*x2*x2*x2;
    x1*x1*x2*x2*x2*x2;
    x1*x2*x2*x2*x2*x2;
    x2*x2*x2*x2*x2*x2];
end