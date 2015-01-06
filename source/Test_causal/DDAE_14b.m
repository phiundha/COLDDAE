% DDAE TEST PROBLEM 14b
% strangeness-index: 3

% advanced problem!!!

phi=@(t) [sin(t);cos(t);-sin(t);-cos(t)];
tau=@(t)1;
t0=0;
x0=phi(t0);

% exact solution
xe = @(t) phi(t);
xed = @(t) [cos(t);-sin(t);-cos(t);sin(t)];
xet = @(t) phi(t-tau(t));

        E=@(t)[
            0   1   0   0
            0   0   1   0
            0   0   0   1
            0   0   0   0
            ];
        A=@(t)[
            1   0   0   0
            0   1   0   0
            0   0   1   0
            0   0   0   1
            ];
        B=@(t)[
            0   0   0   0
            0   0   0   0
            0   0   0   0
            1   0   0   0
            ];
        f=@(t) E(t)*xed(t)-A(t)*xe(t)-B(t)*xet(t);

n =4;
C = zeros(n,n);
D = eye(n);
r = phi(4*pi);
tN = 4*pi;
N = 100;