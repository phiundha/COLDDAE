% DDAE TEST PROBLEM 04
% strangeness-index: 2

clear all; close all; clc

E=@(t) [
    0   1   0
    0   0   1
    0   0   0
    ];
A=@(t) [
    1   0   0
    0   1   0
    0   0   1
    ];
B=@(t) [
    0   1   0
    0   0   0
    0   0   0
    ];

	phi=@(t)[
    exp(t)
    t
    ones(1,length(t))
    ];

tau=@(t)1;
t0=0;
% inconsistent initial vector
%x0=[-3 0 -1]';
x0=[1 0 1]';


% the exact solution and its derivative
xe = @(t)phi(t);
xed = @(t)[
    exp(t)
    1
    0
    ];

f=@(t)E(t)*xed(t) - A(t)*xe(t) - B(t)*xe(t-tau(t));

