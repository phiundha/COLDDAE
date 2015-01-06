% Testing the code colddae on causal systems
% Problems DDAE_05, DDAE_06, DDAE_14 are of advanced type.
% In particular, DDAE_14 shows which components are "badly solved"
% Problem DDAE_09 is intersting - may try BVP solver
% Problem DDAE_12 is not understandable why the error is not good as
% expected

DDAE_01;

tspan = [0,10];
options.Step = 1e-3;
%options.Iter=100;
options.MaxStrIdx = 3;

sourcefolder = genpath('../');
addpath(sourcefolder);

tic
[t,x,info] = colddae(E,A,B,f,tau,phi,tspan,options);
toc

n = length(x(:,1));
X_e = xe(t);

figure(); clf;

for i=1:n
    subplot(2,n,i);
    plot(t,x(i,:))
    ylim([min(x(i,:))*1.1 max(x(i,:)) * 1.1])
    grid
    legend(strcat('x_',num2str(i)));
     
    subplot(2,n,n+i)
    error_i = x(i,:)-X_e(i,:);
    plot(t,error_i)
    legend(strcat('Abs. error of x_',num2str(i)));
end

rmpath(sourcefolder)
    
