% DDAE TEST PROBLEM 15
% strangeness-index: 1

        E=@(t)[
            1   0
            0   0
            ];
        A=@(t)[
            0   0
            1   -1
            ];
        B=@(t)[
            0   0
            0   -1
            ];

        f=@(t)[
            0
            0
            ];

phi=@(t)[1; (t>=0)];

tau=@(t)1;
t0=0;
x0=phi(t0);

% the exact solution
% xe = @(t)phi(t);
xe = [];
