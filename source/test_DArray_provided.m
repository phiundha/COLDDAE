E0=[
    1 0
    0 0
    ];
A0=[
    0 0
    1 0
    ];
B0=[
    0 1
    0 1
    ];

E=@(t)E0;
A=@(t)A0;
B=@(t)B0;

phi = @(t)[sin(t);cos(t)];
dphi = @(t)[cos(t);-sin(t)];
ddphi = @(t)-phi(t);

tau=@(t) exp(-t);
dtau=@(t) -exp(-t);

f=@(t)E0*dphi(t)-A0*phi(t)-B0*phi(t-tau(t));
df = @(t) E0*ddphi(t)-A0*dphi(t)-B0*dphi(t-tau(t))*(1-dtau(t));

% options with derivative array M*z = P*z(t-tau(t))+g.
options1.DArray{1}=@(t)[
     -A0,E0,zeros(2)
     zeros(2),-A0,E0];
options1.DArray{2}=@(t)[
     B0,zeros(2)
     zeros(2),B0*(1-dtau(t))];
options1.DArray{3}=@(t)[f(t);df(t)];
options1.Shift = 1;
options1.StrIdx = 1;
options1.MaxShift = 1;
options1.MaxStrIdx = 1;

% options without d array
options2.Shift = 1;
options2.StrIdx = 1;
options2.MaxShift = 1;
options2.MaxStrIdx = 1;

tspan=[0,10];

% solve with DArray
[t1,x1,info1]=colddae_noncausal(E,A,B,f,tau,phi,tspan,options1)

% solve with finite differences
[t2,x2,info2]=colddae_noncausal(E,A,B,f,tau,phi,tspan,options2)

semilogy(t1,abs(x1-phi(t1)),t2,abs(x2-phi(t2)),'o-')