% This system is noncausal with
% Strangeness_index: 2
% Shift_index: 1

E=@(t) [
    0   1   0
    0   0   1
    0   0   0
    ];
A=@(t) [
    1   0   0
    0   1   0
    0   0   0
    ];
B=@(t) [
    0   1   0
    0   0   0
    0   0   1
    ];

phi=@(t)[
    exp(t)
    t
    ones(1,length(t))
    ];

tau = @(t)t/2+1;
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
