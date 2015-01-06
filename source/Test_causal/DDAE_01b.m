% DDAE TEST PROBLEM 01b
% strangeness-index: 1
% advanced type / good num. sol.
% May be because of the special structure of E,A, so the discretization
% works.

 clear all; close all; clc
        E=@(t)[
            0   1
            0   0
            ];
        A=@(t)[
            1   0
            0   1
            ];
        B=@(t)[
            1   0
            0   1
            ];
        f=@(t)[
            1-exp(t)-exp(t-1)
            -2*t+1
            ];
        
phi=@(t)[
    exp(t)
    t
    ];
tau=1;
t0=0;
x0=phi(t0);
% the exact solution
xe = @(t)phi(t);

