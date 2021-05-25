%% Assignment 1 SC42100 Networked and Distributed Control Systems
clc, clear

a = 4;
b = 5;
c = 7;

A = [a, b+0.5; 0, -c];
B = [0; 1];


%% Question 1.1: constant sampling interval, no delay

% design continuous controller
p = [-2, -3];
K_bar = place(A,B,p);

% check if A is invertable
if rank(A)== length(A)
    print = 'A is invertible';
    disp(print)
else
    error('A is not invertible')
end

% analyse stability for different values of h
hs = 0.1:0.01:0.5;
spec_rad = zeros(length(hs),1);
i = 1;

for h = hs
    F_cl = expm(A*h) - (expm(A*h)-eye(2))/A*B*K_bar;
    spec_rad(i) = max(abs(eig(F_cl)));
    i = i + 1;
end

figure(1)
plot(hs, spec_rad, 'b')
hold on
plot(hs, ones(length(hs),1),'r')
xlabel('h')
ylabel('\rho(F_{cl}(h))')
title('Spectral radius of the closed-loop system')


%% Question 1.2: constant small delay
% analyse stability for different values of h
hs = 0.1:0.01:0.5;
spec_rad = zeros(length(hs),10);
i = 1;

for h = hs
    j = 1;
    for tau = 0:0.1*h:0.9*h
        Fx = expm(A*h);
        Fu = (expm(A*h) - expm(A*(h-tau)))/A*B;
        G1 = (expm(A*(h-tau)) - eye(2))/A*B;
        Ku = 0; % truly static controller, set to 0.16 for the dynamic tuned controller
        F_cl_ext = [Fx Fu; zeros(1,3)] - [G1; 1]*[K_bar Ku];
        spec_rad(i,j) = max(abs(eig(F_cl_ext)));
        j = j + 1;
    end
    i = i + 1;
end

figure(2)
plot(hs, spec_rad(:,1), 'b')
hold on
for x = 2:10
    plot(hs, spec_rad(:,x))
end
plot(hs, ones(length(hs),1),'r')
xlabel('h')
ylabel('\rho(F_{cl}(h, \tau))')
title('Spectral radius of the closed-loop system')
legend('\tau = 0', '\tau = 0.1*h', '\tau = 0.2*h', '\tau = 0.3*h', '\tau = 0.4*h', '\tau = 0.5*h', '\tau = 0.6*h', '\tau = 0.7*h', '\tau = 0.8*h', '\tau = 0.9*h', 'Location', 'southeast')  


%% alternative plot
hs = 0.1:0.01:0.5;
taus = 0:0.01:0.5;
spec_rad = zeros(length(hs),length(taus));
i = 1;

for h = hs
    j = 1;
    for tau = taus
        if tau + 0.0001 < h
            Fx = expm(A*h);
            Fu = (expm(A*h) - expm(A*(h-tau)))/A*B;
            G1 = (expm(A*(h-tau)) - eye(2))/A*B;
            Ku = 0; % truly static controller, set to 0.16 for the dynamic tuned controller
            F_cl_ext = [Fx Fu; zeros(1,3)] - [G1; 1]*[K_bar Ku];
            spec_rad(i,j) = max(abs(eig(F_cl_ext)));
        end
        j = j + 1;
    end
    i = i + 1;
end

figure(3)
hold on
for i = 1:length(hs)
    for j = 1:length(taus)
        if spec_rad(i,j) == 0 
            plot(hs(i), taus(j), 'color',[1,1,1])
        elseif spec_rad(i,j) < 1
            plot(hs(i), taus(j), 'g.')
        else
            plot(hs(i), taus(j), 'r.')
        end
    end
end
xlabel('h')
ylabel('\tau')
title('Combinations of h and \tau for which the system is stable')



%% Question 1.3: constant, larger delays
% analyse stability for different values of h
hs = 0.01:0.001:0.2;
spec_rad = zeros(length(hs),10);
i = 1;

for h = hs
    j = 1;
    for tau = h:0.1*h:1.9*h
        Fx = expm(A*h);
        Fu2 = (expm(A*h) - expm(A*(2*h-tau)))/A*B;
        Fu1 = (expm(A*(2*h-tau)) -eye(2))/A*B;
        Ku = [0 0]; % truly static controller, set to [0.2 0.01] for dynamic tuned controller
        F_cl_ext2 = [Fx Fu1 Fu2; zeros(1,4); 0 0 1 0] - [0; 0; 1; 0]*[K_bar Ku];
        spec_rad(i,j) = max(abs(eig(F_cl_ext2)));
        j = j + 1;
    end
    i = i + 1;
end

figure(3)
plot(hs, spec_rad(:,1), 'b')
hold on
for x = 2:10
    plot(hs, spec_rad(:,x))
end
plot(hs, ones(length(hs),1),'r')
xlabel('h')
ylabel('\rho(F_{cl}(h,\tau))')
title('Spectral radius of the closed-loop system')
legend('\tau = h', '\tau = 1.1*h', '\tau = 1.2*h', '\tau = 1.3*h', '\tau = 1.4*h', '\tau = 1.5*h', '\tau = 1.6*h', '\tau = 1.7*h', '\tau = 1.8*h', '\tau = 1.9*h', 'Location', 'southeast')  


%% alternative plot
hs = 0.01:0.005:0.2;
taus = 0.01:0.005:0.395;
spec_rad = zeros(length(hs),length(taus));
i = 1;

for h = hs
    j = 1;
    for tau = taus
        if tau + 0.001 >= h && tau + 0.0001 < 2*h
            Fx = expm(A*h);
            Fu2 = (expm(A*h) - expm(A*(2*h-tau)))/A*B;
            Fu1 = (expm(A*(2*h-tau)) -eye(2))/A*B;
            Ku = [0.2 0.01]; % truly static controller, set to [0.2 0.01] for dynamic tuned controller
            F_cl_ext2 = [Fx Fu1 Fu2; zeros(1,4); 0 0 1 0] - [0; 0; 1; 0]*[K_bar Ku];
            spec_rad(i,j) = max(abs(eig(F_cl_ext2)));
        end
        j = j + 1;
    end
    i = i + 1;
end

figure(4)
hold on
for i = 1:length(hs)
    for j = 1:length(taus)
        if spec_rad(i,j) == 0 
            plot(hs(i), taus(j), 'color',[1,1,1])
        elseif spec_rad(i,j) < 1
            plot(hs(i), taus(j), 'g.')
        else
            plot(hs(i), taus(j), 'r.')
        end
    end
end
xlabel('h')
ylabel('\tau')
title('Combinations of h and \tau for which the system is stable')


%% Question 1.4: Time-varying delays
% addpath(genpath('yalmip'))
% p = [-2, -3];
% K_bar = place(A,B,p);
% K = [K_bar 0 0];
K = [7.4 3 0 0];

for h = 0.01:0.01:0.15
    % define all closed loop system matrices
    tau1 = 0.2*h;
    Fx = expm(A*h);
    Fu1 = (expm(A*h) - expm(A*(h-tau1)))/A*B;
    G1 = (expm(A*(h-tau1)) - eye(2))/A*B;
    F_cl1= [Fx Fu1 [0;0]; zeros(1,4); 0 0 1 0] - [G1; 1; 0]*K;

    tau2 = 0.5*h;
    Fu2 = (expm(A*h) - expm(A*(h-tau2)))/A*B;
    G2 = (expm(A*(h-tau2)) - eye(2))/A*B;
    F_cl2= [Fx Fu2 [0;0]; zeros(1,4); 0 0 1 0] - [G2; 1; 0]*K;

    tau3 = h;
    Fu3 = (expm(A*(2*h-tau3)) - eye(2))/A*B;
    F_cl3 = [Fx Fu3 [0;0]; zeros(1,4); 0 0 1 0] - [0; 0; 1; 0]*K;

    tau4 = 1.5*h;
    Fu24 = (expm(A*h) - expm(A*(2*h-tau4)))/A*B;
    Fu14 = (expm(A*(2*h-tau4)) - eye(2))/A*B;
    F_cl4 = [Fx Fu14 Fu24; zeros(1,4); 0 0 1 0] - [0; 0; 1; 0]*K;
    
    % optimization
    Q = eye(4);
    P = sdpvar(4,4);
    constraints = [P >= 0, F_cl1'*P*F_cl1 - P + Q <= 0, F_cl2'*P*F_cl2 - P + Q <= 0 ...
         F_cl3'*P*F_cl3 - P + Q <= 0, F_cl4'*P*F_cl4 - P + Q <= 0];
    obj = 0;
    opts = sdpsettings('solver','sedumi','verbose',0);
    res = optimize(constraints, obj, opts);
    if res.problem == 1
        disp('The maximum sampling time h is')
        value(h)
        break
    end
end


%% Question 1.4: Simplified LMI with medium access schedule
% last entry of K will not be used
p = [-3, -4];
K_bar = place(A,B,p);
K = [K_bar 0.2 0];

for h = 0.01:0.01:1
    tau1 = 0.2*h;
    Fx = expm(A*h);
    Fu1 = (expm(A*h) - expm(A*(h-tau1)))/A*B;
    G1 = (expm(A*(h-tau1)) - eye(2))/A*B;
    F_cl1= [Fx Fu1 [0;0]; zeros(1,4); 0 0 1 0] - [G1; 1; 0]*K;

    tau2 = 0.5*h;
    Fu2 = (expm(A*h) - expm(A*(h-tau2)))/A*B;
    G2 = (expm(A*(h-tau2)) - eye(2))/A*B;
    F_cl2= [Fx Fu2 [0;0]; zeros(1,4); 0 0 1 0] - [G2; 1; 0]*K;

    tau3 = h;
    Fu3 = (expm(A*(2*h-tau3)) - eye(2))/A*B;
    F_cl3 = [Fx Fu3 [0;0]; zeros(1,4); 0 0 1 0] - [0; 0; 1; 0]*K;
    
    F_cl_total = F_cl2 * F_cl3 * F_cl1;
    
    Q2 = eye(3);
    P2 = sdpvar(3,3);
    constraints2 = [P2 >= 0, F_cl_total(1:3,1:3)'*P2*F_cl_total(1:3,1:3) - P2 + Q2 <= 0];
    obj = 0;
    opts = sdpsettings('solver','sedumi','verbose',0);
    res2 = optimize(constraints2, obj, opts);
    if res2.problem == 1
        disp('The maximum sampling time h is')
        value(h)
        break
    end
end


