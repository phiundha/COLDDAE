function [t,x,info] = colddae_causal(E,A,B,f,tau,phi,tspan,options)
%COLDDAE_CAUSAL numerical solver for causal linear delay
% differential-algebraic equations of the form
%   E(t)\dot{x}(t)=A(t)x(t)+sum_i B_i(t)x(t-tau_i(t))+f(t)  for t\in(t0,tf]
%             x(t)=phi(t)                                   for t<=t0
% with multiple delay functions tau_i(t)>=tau_min>0 and history function
% phi.
%
% The corresponding DAE Ex = Ax can have strangeness index bigger than
% zero (differentiation index bigger than one). However, note that bigger
% errors might occur if the strangeness index is too big, because
% derivatives are approximated by finite differences.
% The time stepping is based on the method of steps and using Radau IIA
% collocation using constant step size.
%
% @parameters:
%   E,A,B       Coefficients of the DDAE, m-by-n matrix functions.
%   f           m-by-1 vector function.
%   tau         Variable lag, scalar function.
%   phi         n-by-1 history function.
%   tspan       Considered time interval [t0,tf].
%   options     Struct for optional parameters, set by
%               'options.FieldName = FieldValue', see below
%
% @options
%   Iter        The number of time steps, default: 100.
%   Step        The (constant) step size of the Runge-Kutta method, must
%               be smaller than the smallest delay, default:
%               diff(tspan)/100.
%
%   AbsTol      Absolute tolerance, default: 1e-5.
%   RelTol      Relative tolerance, default: 1e-5.
%
%   StrIdx      Lower bound for the strangeness index.
%   MaxStrIdx   Upper bound for the strangeness index.
%
%   InitVal     Initial value, not necessarily consistent.
%   IsConst     A boolean, true if E and A are constant (then the
%               strangeness-free form is computed only once, i.e. the
%               solver needs less computation time), default: false.
%
% @supporting functions:
%   getRegularizedSystem
%   inflateEA
%   inflateB
%   inflatef
%   evalAllHistFuncs
%   nevilleAitken
%
% @return values:
%   t           The discretization of tspan by Iter+1 equidistant points.
%   x           numerical solution at the time nodes in t.
%   info        Struct with information.
%
% @author:
%       Vinh Tho Ma, TU Berlin, mavinh@math.tu-berlin.de
%       Phi Ha, TU Berlin, ha@math.tu-berlin.de

% set the options
%-------------------------------------------------------------------------%
% set missing fields in options
%-------------------------------------------------------------------------%
if ~exist('options','var'),options = {}; end

% bounds for the itererations in the main loop
if ~isfield(options,'MaxIter')   options.MaxIter = 10000; end
if ~isfield(options,'MaxReject') options.MaxReject = 100; end
if ~isfield(options,'MaxCorrect')options.MaxCorrect = 10; end

% initial step size and bounds for the step size
if ~isfield(options,'InitStep')  options.InitStep = diff(tspan)/100; end
if ~isfield(options,'MinStep')   options.MinStep = 0; end
if ~isfield(options,'MaxStep')   options.MaxStep = inf; end

% tolerances
if ~isfield(options,'AbsTol')    options.AbsTol = 1e-5; end
if ~isfield(options,'RelTol')    options.RelTol = 1e-5; end

% the DDAE's indeces (guesses, if unknown)
if ~isfield(options,'StrIdx')    options.StrIdx = 0; end
if ~isfield(options,'MaxStrIdx') options.MaxStrIdx = 3; end

% initial value (not necessarily consistent)
if ~isfield(options,'InitVal')   options.InitVal = phi(tspan(1)); end

% Are E and A constant matrices?
if ~isfield(options,'IsConst')   options.IsConst = false; end

%-------------------------------------------------------------------------%
% defining some parameters
%-------------------------------------------------------------------------%
t0 = tspan(1);
if not(isa(E,'function_handle'))
    error('E must be a function handle.');
end
[m,n] = size(E(0));
h = options.InitStep;
N = options.MaxIter;
x0 = options.InitVal;

% predefining info's fields
info.Solver = 'colddae_causal';
info.Strangeness_index = -1;
info.Advanced = 0;
info.Number_of_differential_eqs = -1;
info.Number_of_algebraic_eqs = -1;
info.Rejected_steps = 0;
info.Computation_time = -1;

%-------------------------------------------------------------------------%
% some more input checks
%-------------------------------------------------------------------------%
% checking tau
if isa(tau,'double')
    tau = @(t)tau;
end
if not(isa(tau,'function_handle'))
    error('Delay tau must be a function handle.');
end
l = numel(tau(0));
% checking A
if not(isa(A,'function_handle'))
    error('A must be a function handle.');
end
if or(size(A(0),1)~=m,size(A(0),2)~=n)
    error('A has wrong size.')
end
% checking B
if not(isa(B,'function_handle'))
    error('B must be a function handle.');
end
if or(size(B(0),1)~=m,size(B(0),2)~=l*n)
    error('B has wrong size.')
end
% checking f
if not(isa(f,'function_handle'))
    error('f must be a function handle.');
end
if or(size(f(0),1)~=m,size(f(0),2)~=1)
    error('f has wrong size.')
end
% checking phi
if or(size(phi(0),1)~=n,size(phi(0),2)~=1)
    error('phi has wrong size.')
end
if options.MaxStrIdx<options.StrIdx
    error('MaxStrIdx must not be less than StrIdx.')
end

% the container for the approximate solution of the DDAE and its stage values
x=nan(3*n,N);
t = nan(1,N);
t(1) = tspan(1);

% find nearest consistent initial value to x0 by determining local
% strangeness-free form at t=t0 and replace x0 by it
hist_funcs = cell(l,1);
for k=1:l
    hist_funcs{k} = phi;
end
BXTAUF = @(s) B(s)*evalAllHistFuncs(hist_funcs,s,tau(s),n,l)+f(s);
[E1,A1,~,A2,g2,mu,Z1,Z2] = getRegularizedSystem(E,A,BXTAUF,t0,options);
if options.IsConst
    options.RegularSystem = {E1,A1,A2,mu,Z1,Z2};
end
info.Strangeness_index = mu;
info.Number_of_differential_eqs = size(E1,1);
info.Number_of_algebraic_eqs = size(A2,1);

% compute nearest consistent initial vector to x0
x(2*n+1:3*n,1)=x0-pinv(A2)*(A2*x0+g2);

%-------------------------------------------------------------------------%
% Time integration.
%-------------------------------------------------------------------------%
tic
for i=1:N-1
    
    absErr = inf;
    relErr = inf;
    
    for j=1:options.MaxReject
        % Estimate the local error by performing a full step and two half
        % steps.
        [x_full,info] = timeStep(E,A,B,f,tau,phi,t(1:i),x(:,1:i),h,options,m,n,info);
        x_half = timeStep(E,A,B,f,tau,phi,t(1:i),x(:,1:i),h/2,options,m,n,info);
        x_half = timeStep(E,A,B,f,tau,phi,[t(1:i),t(i)+h/2],[x(:,1:i),x_half],h/2,options,m,n,info);
        absErr = norm(x_half(2*n+1:3*n)-x_full(2*n+1:3*n));
        relErr = absErr/norm(x_half(2*n+1:3*n));
        % If the error fulfills the prescribed tolerances or the step size
        % is already equal to the minimal step size, then the step is
        % accepted. If not, the step is "rejected", the step size halved,
        % and we repeat the procedure.
        if (absErr<=options.AbsTol && relErr<=options.RelTol) || (h<=options.MinStep)
            info.Rejected_steps = info.Rejected_steps + j-1;
            break;
        end
        h = max(h/2,options.MinStep);
    end
    
    % Use x_half for the approximation at t(i)+c(3)*h = t(i)+h.
    x(:,i+1) = [x_full(1:2*n);x_half(2*n+1:3*n)];
    t(i+1) = t(i) + h;
    
    % Estimate the next step size h.
    h_newAbs = 0.9*h*(options.AbsTol/absErr)^(1/6);
    h_newRel = 0.9*h*(options.RelTol/relErr)^(1/6);
    h = min([h_newAbs,h_newRel,2*h]);
    
    % Impose lower and upper bounds on the step size h.
    h = max(h,options.MinStep);
    h = min([h,options.MaxStep,tspan(2)-t(i+1)]);
    
    if t(i+1)>=tspan(2)
        break
    end
end
% "cutting out" the approximate solution
x=x(2*n+1:3*n,1:i+1);
t=t(1:i+1);
info.Computation_time = toc;

%-------------------------------------------------------------------------%
% Supporting functions.
%-------------------------------------------------------------------------%

function [x_next,info] = timeStep(E,A,B,f,tau,phi,t,x,h,options,m,n,info)

% the data for the RADAU IIA collocation, V is the inverse of A in the
% Butcher tableau
c=[(4-sqrt(6))/10; (4+sqrt(6))/10; 1];
V=[ 3.224744871391589   1.167840084690405  -0.253197264742181
    -3.567840084690405   0.775255128608412   1.053197264742181
    5.531972647421811  -7.531972647421810   5.000000000000000 ];
v0=-V*ones(3,1);

% other parameters
l = numel(tau(t(1)));
tolR = options.RelTol;

% containers for the matrices of the local strangeness-free formulations at
% t_ij = T(i)+h*c(j)
Etij=nan(n,n,3);
Atij=nan(n,n,3);
ftij=nan(n,1,3);

% For each collocation point...
for j=1:3
    TAU = tau(t(end)+c(j)*h);
    hist_funcs = cell(1,l);
    for k = 1:l
        % determine "x(t-tau)" by using the histpry function phi or by
        % using interpolation
        if t(end)+c(j)*h-TAU(k)<t(1)
            hist_funcs{k} = phi;
        else
            % find the biggest time node smaller than t_i+c_j*h-tau
            t_tau_index=find(t(end)+c(j)*h-TAU(k)<t,1)-1;
            if isempty(t_tau_index)
                warning('ONE DELAY BECAME SMALLER THAN THE STEP SIZE. LONG STEPS NOT IMPLEMENTED YET IN solve_''causal_ddae.m''. TERMINATING SOLVING PROCESS.')
                if TAU(k)<options.MinStep
                    error('PLEASE DECREASE THE LOWER BOUND ON THE STEP SIZE.')
                end
                x_next=inf(3*n,1);
                return
            end
            t_tau=t(t_tau_index);
            % if t_i+c_j*h-tau is not a node point, i.e. not in t, then we
            % have to interpolate
            % prepare some data for the interpolation
            % we use a polynomial of degree 3, so we need 4 data points
            x0_tau=x(2*n+1:3*n,t_tau_index);
            X_tau=reshape(x(:,t_tau_index+1),n,3);
            h_tau = t(t_tau_index+1)-t(t_tau_index);
            hist_funcs{k} = @(s) nevilleAitken(t_tau+[0;c]*h_tau,[x0_tau,X_tau],s);
        end
    end

    BXTAUF = @(s) B(s)*evalAllHistFuncs(hist_funcs,s,tau(s),n,l)+f(s);

    % calculate locally regularized (i.e. strangeness-free) form at t =
    % t(i)-c(j)*h
    % we already have E1 and A1, if isConst is TRUE
    if options.IsConst
        E1 = options.RegularSystem{1};
        A1 = options.RegularSystem{2};
        A2 = options.RegularSystem{3};
        mu = options.RegularSystem{4};
        Z1 = options.RegularSystem{5};
        Z2 = options.RegularSystem{6};
        g = zeros((mu+1)*m,1);
        for k = 0:mu
            g(k*m+1:(k+1)*m) = matrixDifferential(BXTAUF,t(end)+c(j)*h,k,tolR,m,1);
        end
        g1 = Z1'*(BXTAUF(t(end)+c(j)*h));
        if numel(Z2)>0
            g2 = Z2'*g;
        else
            g2 = zeros(0,1);
        end
    else
        [E1,A1,g1,A2,g2,mu,~,Z2] = getRegularizedSystem(E,A,BXTAUF,t(end)+c(j)*h,options);
        if mu > info.Strangeness_index
            info.Strangeness_index = mu;
        end
    end

    % check if the derivatives of x(t-tau) vanish in the differential
    % part, if not already known
    if info.Advanced == false
        P = inflateB(B,t(end)+c(j)*h,mu,tolR,m,l*n);
        if numel(Z2)>0
            B_2 = Z2'*P;
            if mu>0
                if (max(max(abs(B_2(:,l*n+1:end))))>tolR*max(max(B_2),1))
                    warning('ACCORDING TO THE CHOSEN TOLERANCE, THE SYSTEM IS VERY LIKELEY OF ADVANCED TYPE, USING THE METHOD OF STEPS MIGHT PRODUCE LARGE ERRORS.')
                    info.Advanced = true;
                end
            end
        end
    end
    Etij(:,:,j)=[E1;zeros(size(A2))];
    Atij(:,:,j)=[A1;A2];
    ftij(:,j)=[g1;g2];
end

% Solve the linear system.
AA=zeros(3*n);
bb=zeros(3*n,1);

for j=1:3
    for k=1:3
        AA((j-1)*n+1:j*n,(k-1)*n+1:k*n)=Etij(:,:,j)/h*V(j,k)-(j==k)*Atij(:,:,j);
    end
    bb((j-1)*n+1:j*n)=ftij(:,j)-Etij(:,:,j)/h*v0(j)*x(2*n+1:3*n,end);
end

% The solution is a vector with length 3*n, it consists of the 3 values
% of the polynomial at the collocation points T(i)+c(j)*h, j=1..3
x_next=AA\bb;


function XTAU = evalAllHistFuncs(hist_funcs,s,TAU,n,l)
% Supporting function for evaluating the functions x(t-tau_i(t)) for i=1:l,
% which are stored as dense output functions in the struct hist_funcs.
XTAU = zeros(l*n,1);
for k = 1:l
    XTAU((k-1)*n+1:k*n,1) = feval(hist_funcs{k},s-TAU(k));
end

function px = nevilleAitken(X,F,x)
% Supporting function for interpolating using the Neville-Aitken algorithm.
n=length(X);
% At first px is a container for the values in the Newton scheme.
px=F;
% Beginning the Newton scheme, see Numerische Mathematik 1.
for i=1:n-1
    for j=1:n-i
        px(:,j)=((x-X(j))*px(:,j+1)-(x-X(j+i))*px(:,j))/(X(j+i)-X(j));
    end
end
px=px(:,1);
function [E_1,A_1,f_1,A_2,f_2,mu,Z1,Z2] = getRegularizedSystem(E,A,f,ti,options)
% A subroutine for regularizing the DAE
%     E(t)\dot{x}(t) = A(t)x(t) f(t),     t\in(t0,tf],
%               x(t) = phi(t),            t<=t0,
% locally at t=ti. The index reduction procedure was taken from:
% -------------------------------------------------------------------------
% P. Kunkel, V. Mehrmann: Differential-Algebraic Equations, Chapter 6.1.
% -------------------------------------------------------------------------

%-------------------------------------------------------------------------%
% defining some parameters
%-------------------------------------------------------------------------%
muMax=options.MaxStrIdx;
tolR=options.RelTol;
mu0=options.StrIdx;
muMax=max(muMax,mu0);
E0 = E(ti);
[m,n]=size(E0);

%-------------------------------------------------------------------------%
% Main loop: Increase mu until we get a regular system.
%-------------------------------------------------------------------------%
if isfield(options,'DArray')
    NM_provided = feval(options.DArray,ti);
end
for mu = mu0:muMax
    % Build the derivative array.
    if isfield(options,'DArray') && mu<=size(NM_provided,1)/m-1
        NM = NM_provided(1:(mu+1)*m,1:(mu+2)*n);
    else
        NM = inflateEA(E,A,ti,mu,tolR);
    end
    M = NM(:,(n+1):end);
    N = -NM(:,1:n);
    % Extract the algebraic equations.
    Z2 = null2(M',tolR);
    A_2 = Z2'*N;
    T2 = null2(A_2,tolR);
    % Check if the number of (linearly independent) algebraic equations
    % a and differential equations d is equal to the number of
    % variables n, if not then continue by increasing mu.
    a = rank(A_2,tolR);
    d = rank(E0*T2,tolR);
    if a+d>=n
        break;
    end
    if mu >= muMax
        error('MAXIMAL NUMBER OF SHIFTS AND STRANGENESS REDUCTION STEPS REACHED. REGULARIZATION OF THE DDAE FAILED.')
    end
end
%-------------------------------------------------------------------------%
% Invariants of the regularized system found, now extracting equations.
%-------------------------------------------------------------------------%
% Remove redundant algebraic equations.
if a>0
    Y2 = orth2(A_2,tolR);
    A_2 = Y2'*A_2;
    % Update the selector of alg. eqs.
    Z2 = Z2*Y2;
end
% Select appropriate differential equations.
Z1 = orth2(E0*T2,tolR);
E_1 = Z1'*E0;
% Extract the algebraic and differential parts for f and the
% differential parts for E, A and B.
g = inflatef(f,ti,mu,tolR,m);
if a>0
    f_2 = Z2'*g;
else
    f_2 = zeros(0,1);
end
A_1 = Z1'*A(ti);
f_1 = Z1'*f(ti);
function NM = inflateEA( E,A,t,mu,tolR )
% Computes the derivative array of (E,A) by differentiating mu times.
%   INPUT
%   -----
%   E       fcn_handle      m-by-n leading matrix function
%   A       fcn_handle      m-by-n matrix function
%   t       double          the time
%   mu      double          the strangeness index
%   tolR    double          the relative tolerance
%
%   OUTPUT
%   ------
%   NM = [-N,M], where
%       M       double((mu+1)*m,(mu+1)*n)   M like in [1]
%       N       double((mu+1)*m,n)          first n columns of N in [1]
%
%   References:
%       [1] P.Kunkel, V. Mehrmann: Differential-Algebraic Equations,
%           chapters 3.1 and 3.2
E0 = E(t);
A0 = A(t);
[m,n] = size(E0);
dE = zeros((mu+1)*m,n);
dA = zeros((mu+1)*m,n);
NM = zeros((mu+1)*m,(mu+2)*n);
dE(1:m,1:n) = E0;
dA(1:m,1:n) = A0;
NM(1:m,n+1:2*n) = E0;
for l = 1:mu
    % Make dE and dA contain all derivatives up to order l.
    dE(l*m+1:(l+1)*m,1:n) = matrixDifferential( E,t,l,tolR,m,n);
    dA(l*m+1:(l+1)*m,1:n) = matrixDifferential( A,t,l,tolR,m,n);
    %Expand M_(l-1) to M_l.
    for j = 0:l-1
        k = l-j;
        NM(l*m+1:(l+1)*m,(j+1)*n+1:(j+2)*n) = nchoosek(l,j)*dE(k*m+1:(k+1)*m,:)-nchoosek(l,j+1) * dA((k-1)*m+1:k*m,:);
    end
    NM(l*m+1:(l+1)*m,(l+1)*n+1:(l+2)*n) = dE(1:m,1:n);
end
NM(:,1:n) = -dA;
function g = inflatef(f,ti,mu,tol,m)
% Builds the vector
%    _        _
%   |          |
%   |   f(ti)  |
%   |   .      |
%   |   f(ti)  |
%   |   ..     |
%   |   f(ti)  |
%   |  ...     |
%   |   f(ti)  |
%   |   .      |
%   |   .      |
%   |   .      |
%   |    (mu)  |
%   |   f (ti) |
%   |_        _|.
%
g = zeros((mu+1)*m,1);
for i = 1:(mu+1)
    g(((i-1)*m+1:((i-1)+1)*m)) = matrixDifferential(f,ti,i-1,tol,m,1);
end
function P = inflateB(B,ti,mu,tol,m,n)
% Builds the matrix
%    _                                         _
%   |                                           |
%   |   B                                       |
%   |   .                                       |
%   |   B    B                                  |
%   |   ..   .                                  |
%   |   B   2B    B                             |
%   |  ...   ..   .                             |
%   |   B   3B   3B    B                        |
%   |   .    .    .    .    .                   |
%   |   .    .    .    .    .    .              |
%   |   .    .    .    .    .    .    .         |
%   |    (mu)                                   |
%   |   B    .    .    .    .    .    .    B    |
%   |_                                         _|.
%
P = zeros((mu+1)*m,(mu+1)*n);
for i = 1:(mu+1)
    P((i-1)*m+1:((i-1)+1)*m,1:n) = matrixDifferential( B,ti,i-1,tol,m,n);
    for j=1:i-1
        k = i-j-1;
        P((i-1)*m+1:((i-1)+1)*m,j*n+1:(j+1)*n) = round(prod((i-j:i-1)./(1:j)))*P(k*m+1:(k+1)*m,1:n);
    end
end
function dA = matrixDifferential(A,t,k,tol,m,n)
% Approximates the time derivative of the (matrix) function A.
eps=0.01;
j=0;
delta=sqrt(eps*max(0.01,abs(t)));
temp=zeros(m,n,k+1);
alpha=tol+1;
while j<2 && alpha>tol
    delta=delta/2;
    dA_old=A(0);
    for i=0:k
        % temp(:,:,i+1)=(-1)^i*nchoosek(k,i)*A(t+(k/2-i)*delta);
        temp(:,:,i+1)=(-1)^i*round(prod(((k-i+1):k)./(1:i)))*A(t+(k/2-i)*delta);
    end
    dA=sum(temp,3)/delta^k;
    alpha=norm(dA-dA_old);
    j=j+1;
end
if min(min(isfinite(dA)))==0
    warning('ERROR IN matrixDifferential.m!')
end
function Z = null2(A,tol)
% Slight modification of MATLAB's null function.
[m,n] = size(A);
[~,S,V]=svd(A,0);
if m > 1
    s = diag(S);
elseif m == 1
    s = S(1);
else s = 0;
end
r = sum(s > max(m,n) * max(s(1),1) * tol);
Z = V(:,r+1:n);
function Q = orth2(A,tol)
% Slight modification of MATLAB's orth function.
if isempty(A)
    Q=A;
    return;
end
[U,S] = svd(A,0);
[m,n] = size(A);
if m > 1, s = diag(S);
elseif m == 1, s = S(1);
else s = 0;
end
r = sum(s > max(m,n) * max(s(1),1) * tol);
Q = U(:,1:r);