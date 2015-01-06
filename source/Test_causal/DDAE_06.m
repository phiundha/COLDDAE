% DDAE TEST PROBLEM 06
% strangeness-index: 2
% System is of advanced type

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
    1   1   1
    1   1   1
    1   1   1
    ];
f=@(t)[
    2*pi*cos(2*pi*t)-sin(2*pi*t)-cos(2*pi*t)
    -2*pi*sin(2*pi*t)-2*sin(2*pi*t)-cos(2*pi*t)
    -2*cos(2*pi*t)-sin(2*pi*t)
    ];

phi=@(t)[
    zeros(1,length(t))
    sin(2*pi*t)
    cos(2*pi*t)
    ];
tau=@(t)1;
t0=0;
x0=[0 0 1]';


% EXACT SOLUTION
xe=@(t)[zeros(1,length(t));sin(2*pi*t);cos(2*pi*t)];