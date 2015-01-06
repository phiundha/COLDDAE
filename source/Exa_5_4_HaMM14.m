% DDAE TEST PROBLEM WITH ADDITIONAL DIFFERENTIAL EQUATIONS
% strangeness-index: 1
% shift index = 1

clear;
E=@(t) [
    1   0   
    0   0   
    ];
A=@(t) [
    0   0  
    1   0  
    ];
B=@(t) [
    0   1   
    0   1   
    ];

phi=@(t)[
    cos(t)
    sin(t)
    ];

tau=@(t) 1-sin(t)/2;
t0=0;
% consistent initial vector
x0=[1 0]';

% the exact solution
xe = phi;
xed = @(t)[-sin(t); cos(t)];

E0=E(0);
A0=A(0);
B0=B(0);
% f=@(t) E0*[ 1; cos(t)] - A0*[t;sin(t)] - B0*[t-tau;sin(t-tau)] ;
f=@(t) E0*xed(t) - A0*xe(t) - B0*xe(t-tau(t)) ;

options.Shift = 1;
options.StrIdx = 1;

options.MaxShift = 1;
options.MaxStrIdx = 1;

options.IsCausal = 0;

tic
[t,x,info]=colddae(E,A,B,f,tau,phi,[0,10],options);
% [t,x]=solve_shifted_lddae(E,A,B,f,tau,phi,[0,10],options);
toc

semilogy(t,abs(x-xe(t))./abs(xe(t)),'o-');
grid
legend('rel. error of x_1(t)','rel. error of x_2(t)');