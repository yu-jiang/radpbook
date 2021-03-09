function b = cost_base(x)
x1 = x(1);
x2 = x(2);
x3 = x(3);

b = [x1*x1;
     2*x1*x2;
     2*x1*x3;
     x2*x2;
     x2*x3;
     x3*x3];
end