% Test example of the paper HaM14
% Example 4.4

E=@(t)[1 0; 0 0];
A=@(t)[0 0; 1 0];
B=@(t)[0 1; 0 0];

%phi = @(t) [sin(t);cos(t)];
%dphi = @(t)[cos(t);-sin(t)];
phi = @(t) [t; exp(t/10)];
dphi = @(t)[1;0.1 * exp(t/10)];
    
%tau = @(t) exp(-t);
tau = @(t)1;

xe = phi;

f=@(t)[1-exp(0.1 * (t-tau(t)));-t];
