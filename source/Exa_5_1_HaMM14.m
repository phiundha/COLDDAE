% Example 1.4.9, taken from 
% C. H.A. Paul,
% A Test Set of Functional Differential Equations
% Numerical Analysis Report No. 243, February 1994
% http://www.maths.manchester.ac.uk/~chris/reports/rep243.pdf

clear all; close all; clc

m=2;

E=@(t)eye(4);
A=@(t)[
    0 0 1 0
    0 0 0 1
    0 -2*m 0 0
    -2*m 0 0 0
    ];
B=@(t)[
    0 0 0 0
    0 0 0 0
    (1+m^2)*(-1)^m 0 0 0
    0 (1+m^2)*(-1)^m 0 0
    ];

f=@(t)zeros(4,1);
tau=@(t)pi;

phi=@(t)[
    sin(t).*cos(m*t)
    cos(t).*sin(m*t)
    cos(t).*cos(m*t)-m*sin(t).*sin(m*t)
    m*cos(t).*cos(m*t)-sin(t).*sin(m*t)
    ];

% The analytical solution
x_e=@(t)[
    sin(t).*cos(m*t)
    cos(t).*sin(m*t)
    cos(t).*cos(m*t)-m*sin(t).*sin(m*t)
    m*cos(t).*cos(m*t)-sin(t).*sin(m*t)
    ];

tspan=[0,10];

%options.MaxStep = 1e-3;
%options.Iter=100;
options.MaxStrIdx=0;
tic
[t,x,info] = colddae(E,A,B,f,tau,phi,tspan,options);
toc

figure(); clf;
subplot(1,2,1)
plot(t,x)
legend('x_1','x_2','x_3','x_4');
xlabel('t')
ylabel('numerical solution')

subplot(1,2,2)
semilogy(t,abs(x-x_e(t))./abs(x_e(t)),'o-')
grid
xlabel('t')
legend('rel. error of x_1(t)','rel. error of x_2(t)','rel. error of x_3(t)','rel. error of x_4(t)')

%suptitle('Numerical solution and relative error')