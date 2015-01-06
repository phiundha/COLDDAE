% DDAE TEST PROBLEM 05
% strangeness-index: 1
% System is of advanced type


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
            1   0
            ];
        f=@(t)[
            (2*pi-1)*cos(2*pi*t)-2*pi*t*sin(2*pi*t)
            -2*sin(2*pi*t)-t*cos(2*pi*t)
            ];

phi=@(t)[
    sin(2*pi*t)
    cos(2*pi*t)
    ];
tau=@(t)1;
t0=0;
x0=phi(t0);

% EXACT SOLUTION
xe=@(t) [sin(2*pi*t);cos(2*pi*t)];