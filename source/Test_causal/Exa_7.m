clear all; close all; clc

P = @(t) [
1 t 0
0 1 t^2
0 0 1
];

E = [
    0 1 0 
    0 0 1
    0 0 0
    ];
A = [
    1 0 0
    0 1 0
    0 0 1
    ];
B = [
    1 1 1 1 1 1
    0 0 0 0 0 0
    0 0 0 0 0 0
    ];

xe = @(t) [cos(t);sin(t);atan(t)];
xed= @(t) [-sin(t);cos(t);1/(1+t.^2)];

phi = xe;
tau = @(t) [1-0.5*cos(t),1-0.5*sin(t)];

f = @(t) P(t)*(E*xed(t)-A*xe(t)-B*[xe(t-1+0.5*cos(t));xe(t-1+0.5*sin(t))]);
E = @(t) P(t)*E;
A = @(t) P(t)*A;
B = @(t) P(t)*B;

tspan=[0,10];

options.Iter=1000;
options.StrIdx = 2;
options.isConst = 0;

tic
[t,x,info]=colddae(E,A,B,f,tau,phi,[0,10],options);
%[t,x] = solve_ddae({E,A,B,f,tau,phi},tspan,options);
toc

semilogy(t,abs(x-xe(t)))
title('absolute error')