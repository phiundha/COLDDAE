% DDAE TEST PROBLEM 10
% strangeness-index: 0

clear all; close all; clc

        E=@(t)[
            1   0
            0   0
            ];
        A=@(t)[
            0   1
            t   1
            ];
        B=@(t)[
            0   0
            0   0
            ];
        f=@(t)[
            1-exp(t)
            -t.^2 - exp(t)
            ];
        
phi=@(t)[
    t
    exp(t)    
    ];
tau=@(t)1;
t0=0;
x0=phi(t0);
% the exact solution
xe = @(t)phi(t);