function dx = simpleSysWrapper(~,x,w)
% System dynamics part
u = w'*Psi_fun(x);
dx = susp_sys(x,u);       % dx as the first 1-4 states of the wrapper
end