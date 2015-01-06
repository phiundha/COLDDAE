% DDAE TEST PROBLEM 11
% Strangeness index 1
% 2 constant delay

clear all; close all; clc

E=@(t)[
    0   t
    0   0
    ];
A=@(t)[
    1   0
    0   1
    ];
B1=@(t)[
    0   1
    0   0
    ];
B2=@(t)[
    -exp(2)   0
    0   0
    ];

B = @(t)[B1(t), B2(t)];

f=@(t)[
    1
    -t
    ];

phi=@(t)[
    exp(t)
    t
    ];
tau=@(t)[1, 2];
t0=0;
x0=phi(t0);
% the exact solution
xe = @(t)phi(t);

