% DDAE TEST PROBLEM 02
% strangeness-index: 1 for t not equal to zero

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
        f=@(t)[
            -exp(t)+1
            -t
            ];

phi=@(t)[
    exp(t)
    t
    ];
tau=@(t)1;
tspan = [0,10];
t0 = tspan(1);

x0=phi(t0);
% the exact solution
xe = @(t)phi(t);
