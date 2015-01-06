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

B1=@(t) [
    0   1   0
    0   0   0
    0   0   0
    ];

B2=@(t) [
    0   0   1
    0   0   0
    0   0   0
    ];

B=@(t)[B1(t),B2(t)];

phi=@(t)[
    exp(t)
    t
    sin(t)
    ];

f=@(t) E(t)* [exp(t) 1 cos(t)]' - A(t) * phi(t) - B1(t)*phi(t-1) - B2(t)*phi(t-t/2-1);

xe = @(t)phi(t);

tau=@(t)[1, t/2+1];

x0=[1 0 0]';

tspan=[0,10];

%options.Iter=100;
options.StrIdx = 2;
options.IsConst = 1;
options.IsCausal = 1;

[t,x,info]=colddae(E,A,B,f,tau,phi,tspan,options);

figure(); clf;

subplot(1,2,1)
plot(t,x)
legend('x_1','x_2','x_3');
xlabel('t')
ylabel('numerical solution')

subplot(1,2,2)
semilogy(t,abs(x-phi(t))./abs(phi(t)),'o-')
grid
xlabel('t')
legend('rel. error of x_1(t)','rel. error of x_2(t)','rel. error of x_3(t)')