% DDAE TEST PROBLEM 01
% strangeness-index: 1

clear all; close all; clc

        E=@(t)[
            1   t
            0   0
            ];
        A=@(t)[
            0   0
            1   t
            ];
        B=@(t)[
            0   1
            0   0
            ];
        f=@(t)[
            2
            -t.^2-t-1
            ];
 
phi=@(t)[
    t+1
    t
    ];

tau = @(t)1;
tspan = [0,10];

t0 = tspan(1);
x0=phi(t0);
% the exact solution
xe = @(t)phi(t);
