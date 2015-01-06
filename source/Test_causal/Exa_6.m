E = [
    0 1 0 
    0 0 1
    0 0 0
    0 1 1
    ];
A = [
    1 0 0
    0 1 0
    0 0 1
    1 1 1
    ];
B = @(t)[
    1 1 1 1 1 1
    0 0 0 0 0 0
    0 0 0 0 0 0
    1 1 1 1 1 1
    ]*exp(sin(t));

xe = @(t) [cos(t);sin(t);atan(t)];
xed= @(t) [-sin(t);cos(t);1/(1+t.^2)];

phi = xe;
tau = @(t) [1-0.5*cos(t),1-0.5*sin(t)];

f = @(t) E*xed(t)-A*xe(t)-B(t)*[xe(t-1+0.5*cos(t));xe(t-1+0.5*sin(t))];
E = @(t) E;
A = @(t) A;
% B = @(t) B;

tspan=[0,10];

options.Iter=1000;
options.StrIdx = 2;
options.isConst = 1;

tic
[t,x,info]=colddae(E,A,B,f,tau,phi,[0,10],options);
% [t,x]=solve_shifted_lddae(E,A,B,f,tau,phi,[0,10],options);
toc

semilogy(t,abs(x-xe(t)))
title('absolute error')