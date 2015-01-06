% DDAE TEST PROBLEM 14
% strangeness-index: 3
% advanced problem!!!

clear all; close all; clc

E=@(t)[
    0   1   0   0
    0   0   1   0
    0   0   0   1
    0   0   0   0
    ];
A=@(t)[
    1   0   0   0
    0   1   0   0
    0   0   1   0
    0   0   0   1
    ];
B=@(t)[
    0   0   0   0
    0   0   0   0
    0   0   0   0
    0   0   0   1
    ];

phi=@(t) [sin(t);cos(t);-sin(t);-cos(t)];
tau=1;
t0=0;
x0=phi(t0);

% exact solution and its derivative
xe = @(t) phi(t);
xed = @(t) [cos(t);-sin(t);-cos(t);sin(t)];
f=@(t)E(t)*xed(t)-A(t)*xe(t)-B(t)*xe(t-tau);

