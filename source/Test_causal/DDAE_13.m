% DDAE TEST PROBLEM 13
% strangeness-index: 1 for t not equal to zero
% 1 time varying delay

clear all; close all; clc

E=@(t)[
    0   t
    0   0
    ];
A=@(t)[
    1   0
    0   1
    ];
B=@(t)[
    0   1
    0   0
    ];        
        
phi=@(t)[
    exp(t)
    t
    ];

tau=@(t) t/2 + 1;

f=@(t)E(t)* [exp(t); 1] - A(t)*phi(t) - B(t) * phi(t-tau(t));

t0=0;
x0=phi(t0);
% the exact solution
xe = @(t)phi(t);