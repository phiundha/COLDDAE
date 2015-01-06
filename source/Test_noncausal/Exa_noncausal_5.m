% System is noncausal of shift index 1
% strangeness index 2
% not of advanced type

E = [0 1 0
     0 0 1
     0 0 0];

A=[1   0   0
   0   1   0
   0   0   0];
    
B=[1   0   0
   0   0   0
   0   0   1];

xe = @(t) [cos(t);sin(t);atan(t)];
xed= @(t) [-sin(t);cos(t);1/(1+t.^2)];

phi = xe;
tau = @(t) 1-0.5*sin(t);

f = @(t) E*xed(t)-A*xe(t)-B*xe(t-tau(t));
E = @(t) E;
A = @(t) A;
B = @(t) B;

tspan=[0,10];
