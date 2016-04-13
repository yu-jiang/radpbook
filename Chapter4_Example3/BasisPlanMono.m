function y = BasisPlanMono(x1,x2)
y = [x1;
    x2;
    x1*x1;
    x1*x2;
    x2*x2;
    x1*x1*x1;
    x1*x1*x2;
    x1*x2*x2;
    x2*x2*x2];
end