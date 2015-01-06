% DDAE TEST PROBLEM 02b
% strangeness-index: 1 for t not equal to zero
% advanced type / bad num. sol.

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
            1   0
            ];
        f=@(t)[
            -exp(t)+1
            -t-exp(t-1)
            ];
 
phi=@(t)[
    exp(t)
    t
    ];
tau=@(t)1;
t0=0;
x0=phi(t0);
% the exact solution
xe = @(t)phi(t);
