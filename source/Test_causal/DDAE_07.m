% DDAE TEST PROBLEM 07
% strangeness-index: 2

clear all; close all; clc
E=@(t) [
    0   1   0
    0   0   t.^2
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
f=@(t)[
    2 - exp(t) - t
    -t
    -1
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

% the exact solution
xe = @(t)phi(t);