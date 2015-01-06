function testDDEsByPaul2()
% Benchmark using the DDEs by C.H. Paul.
%-------------------------------------------------------------------------%
% Executing all function tests
%-------------------------------------------------------------------------%
testEq127()
testEq149()
testEq1412()
testEq214()
testEq224()

%-------------------------------------------------------------------------%
% Constant-delay scalar DDEs
%-------------------------------------------------------------------------%
function testEq127()
E=@(t)1;
A=@(t)0;
B=@(t)-1;
f=@(t)1;
tau=@(t)t-exp(1-1/t);
phi=@log;
alpha=0.5;
tspan=[alpha,10];
options.InitStep=tau(alpha)/2;
options.MinStep=0.1;
options.MaxIdx=0;
options.MaxShift=0;
[t,x,info] = colddae_noncausal(E,A,B,f,tau,phi,tspan,options);
abs_error=max(abs(x-phi(t)));
fprintf('Equation 1.2.7: Abs error %1.2d\n',abs_error)

function testEq149()
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
phi12=@(t)[
    sin(t).*cos(m*t)
    cos(t).*sin(m*t)
    ];
tspan=[0,5];
options.MaxIdx=0;
options.MaxShift=0;
[t,x,info] = colddae_noncausal(E,A,B,f,tau,phi,tspan,options);
abs_error=max(max(abs(x-phi(t))));
fprintf('Equation 1.4.9: Abs error %1.2d\n',abs_error)
figure
semilogy(t,abs(x(1:2,:)-phi12(t))./abs(phi12(t)),'o-')
grid
legend('error of x_1(t)','error of x_2(t)')
title('colddae\_noncausal: relative error')

function testEq1412()
m=2;
E=@(t)eye(2);
A=@(t)[
    0 1 
    0 0
    ];
B=@(t)[
    0 0
    2 0
    ];
f=@(t)[
    0
    exp(sin(t)).*(cos(t).^2-sin(t))-2*exp(-cos(t));
    ];
tau=@(t)pi/2;
phi=@(t)[
    exp(sin(t))
    cos(t).*exp(sin(t))
    ];
tspan=[0,5];
options.MaxIdx=0;
options.MaxShift=0;
[t,x,info] = colddae_noncausal(E,A,B,f,tau,phi,tspan,options);
abs_error=max(max(abs(x-phi(t))));
fprintf('Equation 1.4.12: Abs error %1.2d\n',abs_error)

function testEq214()
E=@(t)[
    1 -1
    0 0];
A=@(t)[
    1 0
    0 1
    ];
B=@(t)[
    0 0
   -1 0
    ];
f=@(t)zeros(2,1);
tau=@(t)1;
phi=@(t)[1;1];
tspan=[0,4];
% options.MinStep=0.1;
options.MaxIdx=0;
options.MaxShift=0;
[t,x,info] = colddae_noncausal(E,A,B,f,tau,phi,tspan,options);
xe1 = @(t) exp(t);
xe2 = @(t) (t-1).*exp(t-1)+xe1(t);
xe3 = @(t) 0.5*(t.^2-2*t).*exp(t-2)+xe2(t);
xe4 = @(t) 1/6*(t.^3-3*t.^2-3*t+9).*exp(t-3)+xe3(t);
xe=@(t) xe1(t).*(0<=t).*(t<=1) + xe2(t).*(1<t).*(t<=2) + xe3(t).*(2<t).*(t<=3) + xe4(t).*(3<t).*(t<=4);
rel_error=max(abs(x(1,:)-xe(t))./abs(xe(t)));
fprintf('Equation 2.1.4: Rel error %1.2d\n',rel_error)
figure
semilogy(t,abs(x(1,:)-xe(t))./abs(xe(t)),'o-')
grid
xlabel('t')
ylabel('x_1(t)')
title('colddae\_noncausal: relative error')

function testEq224()
%some trouble with this one, still need to look into that
E=@(t)[
    1 -1
    0 0
    ];
A=@(t)[
    0 0
    0 1
    ];
B=@(t)[
    0 0
    -1 0
    ];
f=@(t)zeros(2,1);
tau=@(t)0.5-t;
phi=@(t)[exp(-t.^2);exp(-(2*t-0.5).^2)];
tspan=[1/4,1/2-0.0001]; % not solvable at t=1/2
options.MaxStep=0.01;
options.LagTol = 0;
options.MaxIdx=0;
options.MaxShift=0;
[t,x,info] = colddae_noncausal(E,A,B,f,tau,phi,tspan,options);
x_ref = 0.8788261256272320;
abs_error=abs(x(1,end)-x_ref);
fprintf('Equation 2.2.4: Abs error %1.2d\n',abs_error)