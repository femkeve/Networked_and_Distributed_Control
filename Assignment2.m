%% Assignment 2 SC42100 Networked and Distributed Control Systems
clc, clear
% addpath(genpath('yalmip'))

a = 4;
b = 5;
c = 7;

A = [a, b+0.5; 0, -c];
B = [0; 1];

% design continuous controller
p = [-2, -3];
K_bar = place(A,B,p);


%% Question 3.1: Packet dropouts
h = 0.1;
deltas = 0:20;
h_ls = h + deltas*h;

for i = 1:length(deltas)
    Q = eye(2);
    P = sdpvar(2,2);
    constraints = P >= 0;
    h_ls = h + h*deltas(1:i);
    for h_l = h_ls
        F = expm(A*h_l);
        G = expm(A*(h_l-h))*((expm(A*h)-eye(2))/A*B*K_bar);
        F_cl = F - G;
        constraints = [constraints, F_cl'*P*F_cl - P + Q <= 0];
    end
    obj = 0;
    opts = sdpsettings('solver','sedumi','verbose',0);
    res = optimize(constraints, obj, opts);
    if res.problem == 1
        disp('The maximum nr of packet dropouts is')
        value(deltas(i))
        break
    end
end


%% Question 3.2: Bernoulli process
% calculate upper bound p* such that MSS is guaranteed
h = 0.01;
A1 = expm(A*h);                         % packet loss
A0 = A1 - (expm(A*h)-eye(2))/A*B*K_bar; % no packet loss
ps = 0:0.01:1;
result_ber = zeros(length(ps),1);

for p = ps
    P = sdpvar(2,2);
    constraints3 = [P >= 0, P - (1-p)*A0'*P*A0 - p*A1'*P*A1 - 0.0000001 >= 0];
    obj = 0;
    opts = sdpsettings('solver','sedumi', 'verbose', 0);
    res = optimize(constraints3, obj, opts);
    if res.problem == 1
       disp('The upper bound p* on p is')
       value(p)
       break
    end
end

%% Question 3.3: Verify with numerical simulations
h = 0.1;
A1 = expm(A*h);                         % packet loss
A0 = A1 - (expm(A*h)-eye(2))/A*B*K_bar; % no packet loss
x0 = [1;1];
time = 0:0.01:1;
x = zeros(2,length(time));
p_star = 0.3;

for i = 1:length(time)
    r = rand(1);
    if r <= p_star
        x(:,i) = A1*x0;
    else
        x(:,i) = A0*x0;
    end
    x0 = x(:,i);
end

figure(1)
hold on
plot(time, x(1,:), 'r')
title('Simulation of a Bernoulli process')
legend('p* = 0.1', 'p* = 0.3')
xlabel('time (s)')
ylabel('state')


%% Questions 3.4: Bernoulli with ASS
h = 0.1;
A1 = expm(A*h);                         % packet loss
A0 = A1 - (expm(A*h)-eye(2))/A*B*K_bar; % no packet loss
ps = 0.01:0.005:1;
bernoulli = zeros(length(ps),1);

for i = 1:length(ps)
    lambda = ps(i);
    for p = ps
        L = (0.999*lambda^(p-1))^(1/p);
        P = sdpvar(2,2);
        constraints = [P - 0.0001 >= 0, A0'*P*A0 <= lambda*P, ...
            A1'*P*A1 <= L*P];
        obj = 0;
        opts = sdpsettings('solver','sedumi','verbose', 0);
        res = optimize(constraints, obj, opts);
        res.problem
        if res.problem == 1
           disp('For lambda equals')
           value(lambda)
           disp('the upper bound p* on p is')
           value(p)
           bernoulli(i) = p;
           break
        end
    end
end 

figure(2)
plot(ps,bernoulli,'bo')

%% Question 3.5: Gilbert-Elliot: takes long to run! 
% For verification it is suggested to reduce the range of ps
h = 0.1;
ps = 0:0.01:1;
A1 = expm(A*h);                         % packet loss
A0 = A1 - (expm(A*h)-eye(2))/A*B*K_bar; % no packet loss

result_gil = zeros(length(ps));

i = 1;

for p00 = ps
    j = 1;
    for p11 = ps
        P0 = sdpvar(2,2);
        P1 = sdpvar(2,2);
        constraints = [P0 >= 0, P1 >= 0, P0 - p00*A0'*P0*A0 - (1-p00)*A1'*P1*A1 - 0.0000001 >= 0, ...
            P1 - p11*A1'*P1*A1 - (1-p11)*A0'*P0*A0 - 0.0000001 >= 0];
        obj = 0;
        opts = sdpsettings('solver','sedumi', 'verbose', 0);
        res = optimize(constraints, obj, opts);
        result_gil(i,j) = res.problem;
        j = j + 1;
%         if res.problem == 1
%            disp('The upper bound p* on p is')
%            value(p)
%            break
%         end
    end
    i = i + 1;
end


% plot results
figure(3)
hold on
for i2 = 1:length(ps)
    for j2 = 1:length(ps)
        if result_gil(i2,j2) == 0
            plot(ps(i2), ps(j2), 'g.')
        else
            plot(ps(i2), ps(j2), 'r.')
        end
    end
end
xlabel('p_{00}')
ylabel('p_{11}')
title('Combinations of p_{00} and p_{11} for which MSS is guaranteed')


%% Question 4.1: polytopic approximation with four vertices
hs = 0.1:0.01:0.5;
taus = 0.01:0.01:0.5;

jordan = zeros(length(hs),length(taus));

for i = 1:length(hs)
    h = hs(i);
    
    for j = 1: length(taus)
        tau = taus(j);

        if tau - 0.0000001 <= h
            % min value of alpha, for max allowable tau
            alpha_1_min = exp(4*(h-tau)); 
            alpha_2_min = exp(-7*(h-tau));
            
            % max value of alpha, for tau = 0;
            alpha_1_max = exp(4*(h));
            alpha_2_max = exp(-7*(h));

            % calculate matrices            
            F0 = [expm(A*h), [1/14*exp(-7*h)+1/8*exp(4*h); -1/7*exp(-7*h)]; zeros(1,3)];
            F1 = [zeros(3,2),[-1/8; 0; 0]]; 
            F2 = [zeros(3,2),[-1/14; 1/7; 0]];

            G0 = [-11/56; 1/7; 1];
            G1 = [1/8; 0; 0];
            G2 = [1/14; -1/7; 0];

            K = [K_bar 0];

            % solve LMIs
            P = sdpvar(3,3);
            gamma = 0.001;
            constraints = [P - 0.00001 >= 0, ... 
                ((F0+alpha_1_max*F1+alpha_2_max*F2)-(G0+alpha_1_max*G1+alpha_2_max*G2)*K)'*P*((F0+alpha_1_max*F1+alpha_2_max*F2)-(G0+alpha_1_max*G1+alpha_2_max*G2)*K) - P <= -gamma*P, ... 
                ((F0+alpha_1_max*F1+alpha_2_min*F2)-(G0+alpha_1_max*G1+alpha_2_min*G2)*K)'*P*((F0+alpha_1_max*F1+alpha_2_min*F2)-(G0+alpha_1_max*G1+alpha_2_min*G2)*K) - P <= -gamma*P, ... 
                ((F0+alpha_1_min*F1+alpha_2_max*F2)-(G0+alpha_1_min*G1+alpha_2_max*G2)*K)'*P*((F0+alpha_1_min*F1+alpha_2_max*F2)-(G0+alpha_1_min*G1+alpha_2_max*G2)*K) - P <= -gamma*P, ... 
                ((F0+alpha_1_min*F1+alpha_2_min*F2)-(G0+alpha_1_min*G1+alpha_2_min*G2)*K)'*P*((F0+alpha_1_min*F1+alpha_2_min*F2)-(G0+alpha_1_min*G1+alpha_2_min*G2)*K) - P <= -gamma*P];
            obj = 0;
            opts = sdpsettings('solver','sedumi', 'verbose', 0);
            res = optimize(constraints, obj, opts);
            jordan(i,j) = res.problem;
        end
    end
end

% plot
figure(4)
hold on
for i = 1:length(hs)
    for j = 1:length(taus)
        if taus(j) -0.0000001<= hs(i)
            if jordan(i,j) == 0
                plot(hs(i), taus(j), 'g.')
            else
                plot(hs(i), taus(j), 'r.')
            end
        end
    end
end
xlabel('h (s)')
ylabel('\tau (s)')
title('Stability for combinations of (h,\tau)')


%% Question 4.3: more vertices
% divide the range of tau into two pieces
hs = 0.1:0.005:0.5;
taus = 0.01:0.005:0.5;

jordan = zeros(length(hs),length(taus));

for i = 1:length(hs)
    h = hs(i);
    
    for j = 1: length(taus)
        tau = taus(j);

        if tau <= h
            % min value of alpha, for max allowable tau
            alpha_1_min = exp(4*(h-tau));
            alpha_2_min = exp(-7*(h-tau));
            
            % max value of alpha, for tau = 0;
            alpha_1_max = exp(4*(h));
            alpha_2_max = exp(-7*(h));
            
            % mid values
            alpha_1_mid = exp(4*(h-tau/2));
            alpha_2_mid = exp(-7*(h-tau/2));
           
            % calculate matrices
            F0 = [expm(A*h), [1/14*exp(-7*h)+1/8*exp(4*h); -1/7*exp(-7*h)]; zeros(1,3)];
            F1 = [zeros(3,2),[-1/8; 0; 0]]; 
            F2 = [zeros(3,2),[-1/14; 1/7; 0]];

            G0 = [-11/56; 1/7; 1];
            G1 = [1/8; 0; 0];
            G2 = [1/14; -1/7; 0];
            
            K = [K_bar 0];

            % solve LMIs
            P = sdpvar(3,3);
            gamma = 0.001;
            constraints = [P - 0.000001 >= 0, ... 
                ((F0+alpha_1_max*F1+alpha_2_max*F2)-(G0+alpha_1_max*G1+alpha_2_max*G2)*K)'*P*((F0+alpha_1_max*F1+alpha_2_max*F2)-(G0+alpha_1_max*G1+alpha_2_max*G2)*K) - P <= -gamma*P, ... 
                ((F0+alpha_1_max*F1+alpha_2_mid*F2)-(G0+alpha_1_max*G1+alpha_2_mid*G2)*K)'*P*((F0+alpha_1_max*F1+alpha_2_mid*F2)-(G0+alpha_1_max*G1+alpha_2_mid*G2)*K) - P <= -gamma*P, ... 
                ((F0+alpha_1_mid*F1+alpha_2_max*F2)-(G0+alpha_1_mid*G1+alpha_2_max*G2)*K)'*P*((F0+alpha_1_mid*F1+alpha_2_max*F2)-(G0+alpha_1_mid*G1+alpha_2_max*G2)*K) - P <= -gamma*P, ... 
                ((F0+alpha_1_min*F1+alpha_2_mid*F2)-(G0+alpha_1_min*G1+alpha_2_mid*G2)*K)'*P*((F0+alpha_1_min*F1+alpha_2_mid*F2)-(G0+alpha_1_min*G1+alpha_2_mid*G2)*K) - P <= -gamma*P, ... 
                ((F0+alpha_1_mid*F1+alpha_2_min*F2)-(G0+alpha_1_mid*G1+alpha_2_min*G2)*K)'*P*((F0+alpha_1_mid*F1+alpha_2_min*F2)-(G0+alpha_1_mid*G1+alpha_2_min*G2)*K) - P <= -gamma*P, ... 
                ((F0+alpha_1_min*F1+alpha_2_min*F2)-(G0+alpha_1_min*G1+alpha_2_min*G2)*K)'*P*((F0+alpha_1_min*F1+alpha_2_min*F2)-(G0+alpha_1_min*G1+alpha_2_min*G2)*K) - P <= -gamma*P];
            obj = 0;
            opts = sdpsettings('solver','sedumi', 'verbose', 0);
            res = optimize(constraints, obj, opts);
            jordan(i,j) = res.problem;
        end
    end
end

% plot
figure(5)
hold on
for i = 1:length(hs)
    for j = 1:length(taus)
        if taus(j) <= hs(i)
            if jordan(i,j) == 0
                plot(hs(i), taus(j), 'g.')
            else
                plot(hs(i), taus(j), 'r.')
            end
        end
    end
end
xlabel('h (s)')
ylabel('\tau (s)')
title('Stability for combinations of (h, \tau)')


%% Question 5.2
clear, clc

% change value of sigma in the loadpar file
x0s = [[3; 2], [10;10], [-2;10], [10;-40], [1;10], [4;-6], [0; -9], [-2;0], [-1; 8], [5;4], [3;6], [9;4], [-4;3], [6;-9], [-6;-4], [3;-15], [-14;20], [10;-7], [-5;-16], [14;2]];
result = zeros(length(x0s),1); 
diff = zeros(length(x0s)-1,1);
par = loadpar();

for i = 1:length(x0s)
    x0 = x0s(:,i);
    teout = periodic_control(x0); % this function can be changed to periodic_control
    result(i) = length(teout);
    difference = teout(2:end) - teout(1:end-1);
    diff(i) = mean(difference);
end

% Question 5.2
disp(mean(result))

% Question 5.3
disp(mean(diff))

%% Question 5.4: Periodic controller
function teout = periodic_control(y0)
    h_avg = 0.0679;
    sk = y0;
    par = loadpar();
    V_end = 0.1*y0'*par.P*y0;
    tstart = 0;
    tfinal = 1;
    teout = [];
    i = 1;

    while tstart < tfinal && i < 200 %&& sk'*par.P*sk > V_end
        [t,y,] = ode23(@(t,y)f(t,y,sk),[tstart tstart+h_avg],y0); % integrate for h_avg s
        nt = length(t);
        y0 = y(nt,:)';
        if [y0' sk']*par.trig*[y0; sk] > 0        % check condition
            teout = [teout; t(nt)];
            sk = y0;
        end
        tstart = t(nt);
        i = i + 1;
    end
end


%% Question 5: Event-based control functions
function par = loadpar()
a = 4;
b = 5;
c = 7;

par.A = [a, b+0.5; 0, -c];
par.B = [0; 1];

% design continuous controller
p = [-2, -3];
par.K_bar = place(par.A,par.B,p);


Q = eye(2);
par.P = lyap((par.A-par.B*par.K_bar)',Q);
P = par.P;
sigma = 0.9;

A = par.A;
B = par.B;
K = par.K_bar;

par.trig = [A'*P + P*A + sigma*Q,   -P*B*K;...
            -(B*K)'*P,              zeros(2)];

end


function teout = event_control(y0)
tstart = 0;
tfinal = 1;
sk = y0;
teout = [];
par = loadpar();
V_end = 0.1*y0'*par.P*y0;

for i = 1:200
   % Solve until the first terminal event.
   options = odeset('Events',@(t,y)events(t,y,sk));
   [t,y,te,~,~] = ode23(@(t,y)f(t,y,sk),[tstart tfinal],y0,options);

   nt = length(t);
   teout = [teout; te];       
   
   % Set the new initial conditions: the last y(t) is new y0 and new y(sk)
   y0 = y(nt,:)';
   sk = y(nt,:)';
   
   tstart = t(nt);
   
   % if final time is reached before max nr of events, stop
   if tstart >= tfinal
       break
   end
   % if the Lyapunov decrease is satisfied, stop
   %if sk'*par.P*sk <= V_end
   %    break
   %end
end
end


function dydt = f(t,y,sk)
par = loadpar();
dydt = par.A*y - (par.B*par.K_bar)*sk;
end


function [value,isterminal,direction] = events(t,y,sk)
par = loadpar();
% Locate the time when the triggering condition doesn't hold anymore
value = [y' sk']*par.trig*[y; sk]; % trigering condition
isterminal = 1;   % stop the integration
direction = 1;   % event is only located when event function is increasing
end
