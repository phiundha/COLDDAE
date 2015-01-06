% DDAE TEST PROBLEM 09
% strangeness-index: 2
% becoming advanced if we set B(3,2)=1 (or set ep = 1)

clear all; close all; clc

% parameter that describes B(3,2)
ep = 1;

E=@(t) [
    0   1   0
    0   0   1
    0   0   0
    ];
A=@(t) [
    1   t   0
    0   1   0
    0   0   1
    ];

B=@(t) [
    0   1   0
    0   0   0
    0   ep   0
    ];

phi=@(t)[
    exp(t)
    t
    ones(1,length(t))
    ];

dphi = @(t)[exp(t); 1; 0];

tau=@(t)1;

f = @(t) E(t)* dphi(t) - A(t)*phi(t) - B(t)* phi(t-tau(t));

t0=0;
% inconsistent initial vector
%x0=[-3 0 -1]';
x0=[1 0 1]';

% the exact solution
xe = @(t)phi(t);
