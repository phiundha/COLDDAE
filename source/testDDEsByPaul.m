function testDDEsByPaul()
% Benchmark using the DDEs by C.H. Paul.
%-------------------------------------------------------------------------%
% Executing all function tests
%-------------------------------------------------------------------------%
testEq111()
testEq112()
testEq113()
testEq114()
testEq115()
testEq116()
testEq117()
testEq118()
testEq119()
testEq1110()
testEq1111()
testEq1112()
testEq1113()
% testEq121()
% testEq123()
%-------------------------------------------------------------------------%
% Constant-delay scalar DDEs
%-------------------------------------------------------------------------%
function testEq111()
% linear, constant delay
A=1;B=1;C=1; % parameters as in the paper!
x_exact = @(t) C*sum(A.^(0:(round(t/B)+1)).*(t-((0:(round(t/B)+1))-1)*B).^(0:(round(t/B)+1))./factorial(0:(round(t/B)+1)));
options.isConst=1;
[t,x] = solve_lddae(@(t)1,@(t)0,@(t)A,@(t)0,@(t)B,@(t)C,[0,10],options);
rel_error=abs(14640251/44800-x(end))/14640251*44800;
fprintf('Equation 1.1.1: %1.2d\n',rel_error)
function testEq112()
% linear, constant delay
E=@(t)1;
A=@(t)0;
B=@(t)-1;
f=@(t)0;
tau=@(t)pi/2;
phi=@sin;
tspan=[0,10];
options.isConst=1;
[t,x] = solve_lddae(E,A,B,f,tau,phi,tspan,options);
rel_error=max(abs((x-phi(t))./phi(t)));
fprintf('Equation 1.1.2: %1.2d\n',rel_error)
function testEq113()
% nonlinear, constant delay
lambda=3;
F1=@(t,x,xd,xt)-xd-lambda*xt*(1+x);
F2=@(t,x,xt)[];
tau=@(t)1;
phi=@(t)t;
tspan=[0,20];
options.StepSize=0.1;
[t,x] = solve_nddae(F1,F2,tau,phi,tspan,options);
x20=4.671437497493366;
rel_error=abs(x(end)-x20)/x20;
fprintf('Equation 1.1.3: %1.2d\n',rel_error)
function testEq114()
% nonlinear, constant delay
lambda=1.4;
F1=@(t,x,xd,xt)-xd+(lambda-xt)*x;
F2=@(t,x,xt)[];
tau=@(t)1;
phi=@(t)0.01;
tspan=[0,10];
options.StepSize=0.1;
[t,x] = solve_nddae(F1,F2,tau,phi,tspan,options);
%reference solution
x10=1.367208017754907;
rel_error=abs(x(end)-x10)/x10;
fprintf('Equation 1.1.4: %1.2d\n',rel_error)
function testEq115()
% linear, constant delay
E=@(t)1;
A=@(t)1;
B=@(t)1;
f=@(t)3*cos(t)+5*sin(t);
tau=@(t)pi;
phi=@(t)3*sin(t)-5*cos(t);
tspan=[0,10];
options.isConst=1;
[t,x] = solve_lddae(E,A,B,f,tau,phi,tspan,options);
rel_error=max(abs((x-phi(t))./phi(t)));
fprintf('Equation 1.1.5: %1.2d\n',rel_error)
function testEq116()
% nonlinear, multiple constant delay
F1=@(t,x,xd,xt)-xd-xt(1)+xt(2)-xt(3)*xt(4);
F2=@(t,x,xt)[];
tau=@(t)1:4;
phi=@(t)t<0;
tspan=[0,5];
options={};
[t,x] = solve_nddae(F1,F2,tau,phi,tspan,options);
% exact solution
xe1 = @(t) (0<=t & t<1).*polyval([-1,0],t);
xe2 = @(t) (1<=t & t<2).*polyval([1/2,-1,-1/2],t);
xe3 = @(t) (2<=t & t<3).*polyval([-1/6,1/2,0,-7/6],t);
xe4 = @(t) (3<=t & t<4).*polyval([1/24,-1/6,-1/4,1,-19/24],t);
xe5 = @(t) (4<=t & t<=5).*polyval([-1/120,1/6,-5/3,109/12,-24,2689/120],t);
xe = @(t) xe1(t)+xe2(t)+xe3(t)+xe4(t)+xe5(t);
rel_error=max(abs((x-xe(t))./xe(t)));
fprintf('Equation 1.1.6: %1.2d\n',rel_error)
function testEq117()
% nonlinear, constant delay
alpha   = 1/2;
lambda  = 2;
theta   = 1;
tau     = @(t)2;
gamma   = -1;
m       =7;
F1=@(t,x,xd,xt) -xd+(lambda*theta^m*xt)/(theta^m+xt^m)+gamma*x;
F2=@(t,x,xt)[];
phi=@(t)alpha;
tspan=[0,20];
options.StepSize=0.1;
[t,x] = solve_nddae(F1,F2,tau,phi,tspan,options);
% reference solution
x20 = 1.202617750066138;
rel_error=abs((x(end)-x20)/x20);
fprintf('Equation 1.1.7: %1.2d\n',rel_error)
function testEq118()
% nonlinear, constant delay
alpha   = 999;
K       = 1000;
mu      = 1;
q       = 1.5;
z       = 2.5;
tau     = @(t)2;
F1=@(t,x,xd,xt) -xd-mu*x+mu*xt*(1+q*(1-(xt/K)^z));
F2=@(t,x,xt)[];
phi=@(t)alpha;
tspan=[0,20];
options.StepSize=0.1;
[t,x] = solve_nddae(F1,F2,tau,phi,tspan,options);
% reference solution
x20 = 820.185301965284;
rel_error=abs((x(end)-x20)/x20);
fprintf('Equation 1.1.8: %1.2d\n',rel_error)
function testEq119()
% linear, constant delay
E=@(t)1;
A=@(t)5;
B=@(t)1;
f=@(t)0;
tau=1;
phi=@(t)5;
tspan=[0,2];
options.isConst=1;
[t,x] = solve_lddae(E,A,B,f,tau,phi,tspan,options);
xe=@(t)(t>=0 & t<1).*(6*exp(5*t)-1)+(6*(exp(5)+t-6/5).*exp(5*t-5)+1/5).*(t>=1 & t<=2);
rel_error=max(abs((x-xe(t))./xe(t)));
fprintf('Equation 1.1.9: %1.2d\n',rel_error)
function testEq1110()
% nonlinear, constant delay
tau = @(t)pi;
F1=@(t,x,xd,xt) -xd+xt*x;
F2=@(t,x,xt)[];
phi=@(t)-2*(-pi/2<=t & t<0)-1*(t==0);
tspan=[0,6];
options.StepSize=6/100;
[t,x] = solve_nddae(F1,F2,tau,phi,tspan,options);
% exact solution
xe1 = @(t) -1.*(0<=t & t<pi/2);
xe2 = @(t) -exp(pi-2*t).*(pi/2<=t & t<pi);
xe3 = @(t) -exp(-t).*(pi<=t & t<3*pi/2);
xe4 = @(t) -exp(-3/2*pi+1/2*(exp(3*pi-2*t)-1)).*(3*pi/2<=t & t<=6);
xe  = @(t) xe1(t)+xe2(t)+xe3(t)+xe4(t);
rel_error=max(abs((x-xe(t))./xe(t)));
fprintf('Equation 1.1.10: %1.2d\n',rel_error)
function testEq1111()
% nonlinear, constant delay
tau = @(t)pi;
F1=@(t,x,xd,xt) -xd-xt*x;
F2=@(t,x,xt)[];
phi=@(t)-2*(-pi/2<=t & t<0)-1*(t==0);
tspan=[0,6];
options.StepSize=6/100;
[t,x] = solve_nddae(F1,F2,tau,phi,tspan,options);
% exact solution
xe1 = @(t) -1.*(0<=t & t<pi/2);
xe2 = @(t) -exp(2*t-pi).*(pi/2<=t & t<pi);
xe3 = @(t) -exp(t).*(pi<=t & t<3*pi/2);
xe4 = @(t) -exp(3/2*pi+1/2*(exp(2*t-3*pi)-1)).*(3*pi/2<=t & t<=6);
xe  = @(t) xe1(t)+xe2(t)+xe3(t)+xe4(t);
rel_error=max(abs((x-xe(t))./xe(t)));
fprintf('Equation 1.1.11: %1.2d\n',rel_error)
function testEq1112()
%linear, constant delay
E=@(t)1;
A=@(t)1;
B=@(t)1;
f=@(t)0;
tau=@(t)1;
phi=@(t)-1/3<=t ;
tspan=[0,8/3];
options.isConst=1; % DDE with constant coefficients
[t,x] = solve_lddae(E,A,B,f,tau,phi,tspan,options);
% exact solution
c1 = 1+exp(-2/3);
c2 = c1-2*exp(-1);
c3 = 5/3*exp(-1)*(1-c1)+c2-exp(-5/3);
c4 = exp(-2)+c3+2*(c1-c2)*exp(-1);
xe1 = @(t) exp(t).*(0<=t & t<2/3);
xe2 = @(t) (c1*exp(t)-1).*(2/3<=t & t<1);
xe3 = @(t) (t.*exp(t-1)+c2*exp(t)).*(1<=t & t<5/3);
xe4 = @(t) (1+c1*t.*exp(t-1)+c3*exp(t)).*(5/3<=t & t<2);
xe5 = @(t) ((1/2*t.^2-t).*exp(t-2)+c2*t.*exp(t-1)+c4*exp(t)).*(2<=t & t<=8/3);
xe  = @(t) xe1(t)+xe2(t)+xe3(t)+xe4(t)+xe5(t);
rel_error=max(abs((x-xe(t))./xe(t)));
fprintf('Equation 1.1.12: %1.2d\n',rel_error)
function testEq1113()
% nonlinear, constant delay
r = 3.5;
M = 19;
tau = @(t)37/50;
F1=@(t,x,xd,xt) -xd+r*x*(1-xt/M);
F2=@(t,x,xt)[];
phi=@(t)19.001;
tspan=[0,40];
options.StepSize=0.1;
[t,x] = solve_nddae(F1,F2,tau,phi,tspan,options);
% reference solution
x40 = 24.7645486398692;
rel_error=abs((x(end)-x40)./x40);
fprintf('Equation 1.1.13: %1.2d\n',rel_error)

%-------------------------------------------------------------------------%
% Varying-delay scalar DDEs
%-------------------------------------------------------------------------%
function testEq121()
%linear, time-dependent delay
a=3;
b=19/20;
c=11/10;
d=1;
E=@(t)1;
A=@(t)0;
B=@(t)a;
f=@(t)0;
tau=@(t)t-b*t^c;
phi=@(t)d ;
tspan=[0,b^(1/1-c)];
options.MaxStrIdx=0;
options.MaxShift=0;
[t,x] = coldedaNonCausal(E,A,B,f,tau,phi,tspan,options);
% exact solution
% NOT SOLVABLE WITH OUR CURRENT METHOD OF STEPS, SINCE THE DELAY IS SMALLER
% THAN THE STEP SIZE
x1 = 91.22490537957470909;
rel_error=abs(x(end)-x1)/x1;
fprintf('Equation 1.2.1: %1.2d\n',rel_error)

function testEq123()
%linear, time-dependent delay
E=@(t)1;
A=@(t)0;
B=@(t)1;
f=@(t)0;
tau=@(t)1-t;
phi=@(t)1 ;
tspan=[0,1]; % solve only until 0.9 because the delay becomes smaller as t goes to 1
options.MaxStrIdx=0;
options.MaxShift=0;
[t,x] = coldedaNonCausal(E,A,B,f,tau,phi,tspan,options);
% reference solution
x1 = 2.271492555500812;
rel_error=abs(x(end)-x1)/x1;
fprintf('Equation 1.1.13: %1.2d\n',rel_error)
plot(t,x,'o-')