%% SC42100 Networked and Distributed Control Systems Assignment 4
% Femke van Engen, 4552687
clc, clear

%% Dynamics and initial positions for each plane
% state vector is defined as (x,y,xdot,ydot)

A1 =    [1 0 2 0; 0 1 0 2; 0 0 3 0; 0 0 0 3];
B1 =    [2 0;0 2;3 0;0 3];
x01 =   [-10;10;-1;1];

A2 =    [1 0 3 0; 0 1 0 3; 0 0 7 0; 0 0 0 7];
B2 =    [3 0; 0 3; 7 0; 0 7];
x02 =   [10;10;1;1];

A3 =    [1 0 1 0; 0 1 0 1; 0 0 1.1 0; 0 0 0 1.1];
B3 =    [1 0; 0 1; 1.1 0; 0 1.1];
x03 =   [10;-10;1;-1];

A4 =    [1 0 6 0; 0 1 0 6; 0 0 20 0; 0 0 0 20];
B4 =    [6 0;0 6;20 0; 0 20];
x04 =   [-10;-10;-1;-1];

Tfinal = 5;
umax = 100;

time = 0:Tfinal-1;
N = length(time);


%% Centralized solution

A = blkdiag(A1,A2,A3,A4);
B = blkdiag(B1,B2,B3,B4);
x0 = [x01; x02; x03; x04];

% prediction matrices to find the states from initial conditions (T) and the inputs (S)

dim_x = length(A);
dim_u = size(B,2);

T = zeros(dim_x*(N),dim_x);
for k=1:N
    T((k-1)*dim_x+1:k*dim_x,:)=A^k;
end

S = zeros(dim_x*(N),dim_u*(N));
for k=1:N
    for i=0:k-1
        S((k-1)*dim_x+1:k*dim_x,i*dim_u+1:(i+1)*dim_u)=A^(k-1-i)*B;
    end
end

% solve problem with cvx with the SDPT3 solver
cvx_clear;
cvx_begin
    % optimization variable
    variables u(dim_u*N) xfinal(4)
    % objective function
    minimize( (S*u+T*x0)'*(S*u+T*x0) + u'*u )
    % constraints
    subject to 
         S(end-15:end,:)*u == [xfinal;xfinal;xfinal;xfinal] - T(end-15:end,:)*x0;
         norm(u,Inf) <= umax/Tfinal;
cvx_end


%% Plot trajectories (x or y) for all planes
state = S*u+T*x0;
plt_x = [];
plt_y = [];

for i = 0:N-1
    plt_x = [plt_x; state(i*dim_x+1),state(i*dim_x+5),state(i*dim_x+9),state(i*dim_x+13)];
end
for i = 0:N-1
    plt_y = [plt_y; state(i*dim_x+2),state(i*dim_x+6),state(i*dim_x+10),state(i*dim_x+14)];
end

figure(1)
plot (plt_x)
legend('plane 1', 'plane 2', 'plane 3', 'plane 4')
title('Position of the planes in the x-direction')
xlabel('Time step')
ylabel('Position')

figure(2)
plot (plt_y)
legend('plane 1', 'plane 2', 'plane 3', 'plane 4')
title('Position of the planes in the y-direction')
xlabel('Time step')
ylabel('Position')

figure(3)
plot(plt_x,plt_y,'o-')
legend('plane 1', 'plane 2', 'plane 3', 'plane 4')
title('Position of the planes over time')
xlabel('Position x-direction')
ylabel('Position y-direction')


%% Dual decomposition-based algorithm using the projected subgradient-method

% Compute the prediction matrices for every plane
T1 = zeros(4*(N),4);
for k=1:N
    T1((k-1)*4+1:k*4,:)=A1^k;
end

S1 = zeros(4*(N),2*(N));
for k=1:N
    for i=0:k-1
        S1((k-1)*4+1:k*4,i*2+1:(i+1)*2)=A1^(k-1-i)*B1;
    end
end

T2 = zeros(4*(N),4);
for k=1:N
    T2((k-1)*4+1:k*4,:)=A2^k;
end

S2 = zeros(4*(N),2*(N));
for k=1:N
    for i=0:k-1
        S2((k-1)*4+1:k*4,i*2+1:(i+1)*2)=A2^(k-1-i)*B2;
    end
end

T3 = zeros(4*(N),4);
for k=1:N
    T3((k-1)*4+1:k*4,:)=A3^k;
end

S3 = zeros(4*(N),2*(N));
for k=1:N
    for i=0:k-1
        S3((k-1)*4+1:k*4,i*2+1:(i+1)*2)=A3^(k-1-i)*B3;
    end
end

T4 = zeros(4*(N),4);
for k=1:N
    T4((k-1)*4+1:k*4,:)=A4^k;
end

S4 = zeros(4*(N),2*(N));
for k=1:N
    for i=0:k-1
        S4((k-1)*4+1:k*4,i*2+1:(i+1)*2)=A4^(k-1-i)*B4;
    end
end


%% Dual ascent method

iter = 700;
alphas = 0.1;
%alphas = [0.5, 0.1, 0.05, 0.01, 0.005];  % uncomment to run algorithm for multiple step sizes
errors_total = [];
beta = 1;


for a = 1:length(alphas)
    
    % initialize lambdas
    lambda1 = zeros(4,iter);
    lambda2 = zeros(4,iter);
    lambda3 = zeros(4,iter);
    lambda4 = zeros(4,iter);
    subs1 = zeros(4,iter);
    subs2 = zeros(4,iter);
    subs3 = zeros(4,iter);
    subs4 = zeros(4,iter);
    errors = zeros(4,iter);
    alpha = alphas(a);

    for k = 1:iter
        %alpharates = [1/(10+k); 0.1/sqrt(k)]; % uncomment these lines to run for different step size rates
        %alpha = alpharates(a);
        
        % node 1 determines u1 that minimizes the objective and returns xf1
        cvx_clear
        cvx_begin quiet
            variable u1(2*N) 
            minimize ((S1*u1+T1*x01)'*(S1*u1+T1*x01) + u1'*u1 + (lambda1(:,k)-lambda4(:,k))'*(S1(end-3:end,:)*u1+T1(end-3:end,:)*x01))
            subject to
                norm(u1,Inf) <= umax/Tfinal;
        cvx_end

        xf1 = S1(end-3:end,:)*u1+T1(end-3:end,:)*x01;

        % node 2 
        cvx_clear
        cvx_begin quiet
            variable u2(2*N) 
            minimize ((S2*u2+T2*x02)'*(S2*u2+T2*x02) + u2'*u2 + (lambda2(:,k)-lambda1(:,k))'*(S2(end-3:end,:)*u2+T2(end-3:end,:)*x02) )
            subject to
                norm(u2,Inf) <= umax/Tfinal;
        cvx_end 

        xf2 = S2(end-3:end,:)*u2+T2(end-3:end,:)*x02 ;

        % node 3 
        cvx_clear
        cvx_begin quiet
            variable u3(2*N) 
            minimize ((S3*u3+T3*x03)'*(S3*u3+T3*x03) + u3'*u3 + (lambda3(:,k)- lambda2(:,k))'*(S3(end-3:end,:)*u3+T3(end-3:end,:)*x03) )
            subject to
                norm(u3,Inf) <= umax/Tfinal;
        cvx_end 

        xf3 = S3(end-3:end,:)*u3+T3(end-3:end,:)*x03;

        % node 4
        cvx_clear
        cvx_begin quiet
            variable u4(2*N) 
            minimize ((S4*u4+T4*x04)'*(S4*u4+T4*x04) + u4'*u4 + (lambda4(:,k)- lambda3(:,k))'*(S4(end-3:end,:)*u4+T4(end-3:end,:)*x04) )
            subject to
                norm(u4,Inf) <= umax/Tfinal;
        cvx_end 

        xf4 = S4(end-3:end,:)*u4+T4(end-3:end,:)*x04 ;

        % store value of each subgradient and their errors
        subs1(:,k) = xf2-xf1;
        subs2(:,k) = xf3-xf2;
        subs3(:,k) = xf4-xf3;
        subs4(:,k) = xf1-xf4;
        errors(:,k) = [norm(subs1(:,k),1); norm(subs2(:,k),1); norm(subs3(:,k),1); norm(subs4(:,k),1)];
        
        % subgradient iteration to update all the lambdas, then a projection step 
        % the first option for lambda is regular dual ascent, second is Nesterov
        
        %if k > 2 % uncomment this loop for Nesterov
            lambda1(:,k+1) = lambda1(:,k) - alpha*(xf2-xf1);
            %lambda1(:,k+1) = lambda1(:,k) + beta*(lambda1(:,k)-lambda1(:,k-1)) - alpha*(subs1(:,k)+ beta*(subs1(:,k)-subs1(:,k-1)));
            for i = 1:length(lambda1(:,k+1))
                lambda1(i,k+1) = max(lambda1(i,k+1),0);
            end

            lambda2(:,k+1) = lambda2(:,k) - alpha*(xf3-xf2);
            %lambda2(:,k+1) = lambda2(:,k) + beta*(lambda2(:,k)-lambda2(:,k-1)) - alpha*(subs2(:,k)+ beta*(subs2(:,k)-subs2(:,k-1)));
            for i = 1:length(lambda2(:,k+1))
                lambda2(i,k+1) = max(lambda2(i,k+1),0);
            end

            lambda3(:,k+1) = lambda3(:,k) - alpha*(xf4-xf3);
            %lambda3(:,k+1) = lambda3(:,k) + beta*(lambda3(:,k)-lambda3(:,k-1)) - alpha*(subs3(:,k)+ beta*(subs3(:,k)-subs3(:,k-1)));
            for i = 1:length(lambda3(:,k+1))
                lambda3(i,k+1) = max(lambda3(i,k+1),0);
            end

            lambda4(:,k+1) = lambda4(:,k) - alpha*(xf1-xf4);
            %lambda4(:,k+1) = lambda4(:,k) + beta*(lambda4(:,k)-lambda4(:,k-1)) - alpha*(subs4(:,k)+ beta*(subs4(:,k)-subs4(:,k-1)));
            for i = 1:length(lambda4(:,k+1))
                lambda4(i,k+1) = max(lambda4(i,k+1),0);
            end
        %end
    end
    
    errors_total = [errors_total; errors];
end


%% Combined consensus/ incremental subgradient method

iter = 5000;
alpha = 0.0001;
errors_total = [];
xf1 = zeros(4,iter);
xf2 = zeros(4,iter);
xf3 = zeros(4,iter);
xf4 = zeros(4,iter);
lambdas = zeros(4,iter+1);

% consensus matrix
W0 = [0.75 0.25 0 0; 0.25 0.5 0.25 0; 0 0.25 0.5 0.25; 0 0 0.25 0.75];
phi = [0,1,10,100];

% calculate the max and min feasible xf
xf1_min = S1(end-3:end,:)*(-umax/Tfinal*ones(10,1))+T1(end-3:end,:)*x01;
xf1_max = S1(end-3:end,:)*(umax/Tfinal*ones(10,1))+T1(end-3:end,:)*x01;
xf2_min = S2(end-3:end,:)*(-umax/Tfinal*ones(10,1))+T2(end-3:end,:)*x02;
xf2_max = S2(end-3:end,:)*(umax/Tfinal*ones(10,1))+T2(end-3:end,:)*x02;
xf3_min = S3(end-3:end,:)*(-umax/Tfinal*ones(10,1))+T3(end-3:end,:)*x03;
xf3_max = S3(end-3:end,:)*(umax/Tfinal*ones(10,1))+T3(end-3:end,:)*x03;
xf4_min = S4(end-3:end,:)*(-umax/Tfinal*ones(10,1))+T4(end-3:end,:)*x04;
xf4_max = S4(end-3:end,:)*(umax/Tfinal*ones(10,1))+T4(end-3:end,:)*x04;

% define local objectives and constraints
H1 = 2*(S1'*S1+eye(size(S1'*S1)));
f1 = 2*S1'*T1*x01;
H2 = 2*(S2'*S2+eye(size(S2'*S2)));
f2 = 2*S2'*T2*x02;
H3 = 2*(S3'*S3+eye(size(S3'*S3)));
f3 = 2*S3'*T3*x03;
H4 = 2*(S4'*S4+eye(size(S4'*S4)));
f4 = 2*S4'*T4*x04;

lb = -umax/Tfinal*ones(10,1);
ub = umax/Tfinal*ones(10,1);
Aeq1 = S1(end-3:end,:);
Aeq2 = S2(end-3:end,:);
Aeq3 = S3(end-3:end,:);
Aeq4 = S4(end-3:end,:);


for a = 1:length(phi)
    
    % initialize lambdas
    errors = zeros(1,iter);
    W = W0^phi(a);

    for k = 1:iter        
        % node 1
        beq1 = xf1(:,k)-T1(end-3:end,:)*x01;
        [u1,fval1,~,~,lambda1] = quadprog(H1,f1,[],[],Aeq1,beq1,lb,ub,[],[]);
        x1final = S1(end-3:end,:)*u1+T1(end-3:end,:)*x01;
        
        % node 2 
        beq2 = xf2(:,k)-T2(end-3:end,:)*x02;
        [u2,fval2,~,~,lambda2] = quadprog(H2,f2,[],[],Aeq2,beq2,lb,ub,[],[]);
        x2final = S2(end-3:end,:)*u2+T2(end-3:end,:)*x02;
        
        % node 3 
        beq3 = xf3(:,k)-T3(end-3:end,:)*x03;
        [u3,fval3,~,~,lambda3] = quadprog(H3,f3,[],[],Aeq3,beq3,lb,ub,[],[]);
        x3final = S3(end-3:end,:)*u3+T3(end-3:end,:)*x03;
        
        % node 4
        beq4 = xf4(:,k)-T4(end-3:end,:)*x04;
        [u4,fval4,~,~,lambda4] = quadprog(H4,f4,[],[],Aeq4,beq4,lb,ub,[],[]);
        x4final = S4(end-3:end,:)*u4+T4(end-3:end,:)*x04;
        
        % store value of the errors
        errors(k) = norm((x1final+x2final+x3final+x4final)/4-xfinal,1);
        
        % update xfinal values per node and project on feasible set
        xf1(:,k+1) = W(1,1)*(x1final-alpha*-lambda1.eqlin) + W(1,2)*(x2final-alpha*-lambda2.eqlin) ...
                   + W(1,3)*(x3final-alpha*-lambda3.eqlin) + W(1,4)*(x4final-alpha*-lambda4.eqlin);
        for i = 1:length(xf1(:,k+1))
            xf1(i,k+1) = max(xf1(i,k+1), xf1_min(i));
            xf1(i,k+1) = min(xf1(i,k+1), xf1_max(i));
        end
        
        xf2(:,k+1) = W(2,1)*(x1final-alpha*-lambda1.eqlin) + W(2,2)*(x2final-alpha*-lambda2.eqlin) ...
                   + W(2,3)*(x3final-alpha*-lambda3.eqlin) + W(2,4)*(x4final-alpha*-lambda4.eqlin);
        for i = 1:length(xf2(:,k+1))
            xf2(i,k+1) = max(xf2(i,k+1), xf2_min(i));
            xf2(i,k+1) = min(xf2(i,k+1), xf2_max(i));
        end
        
        xf3(:,k+1) = W(3,1)*(x1final-alpha*-lambda1.eqlin) + W(3,2)*(x2final-alpha*-lambda2.eqlin) ...
                   + W(3,3)*(x3final-alpha*-lambda3.eqlin) + W(3,4)*(x4final-alpha*-lambda4.eqlin);
        for i = 1:length(xf3(:,k+1))
            xf3(i,k+1) = max(xf3(i,k+1), xf3_min(i));
            xf3(i,k+1) = min(xf3(i,k+1), xf3_max(i));
        end
        
        xf4(:,k+1) = W(4,1)*(x1final-alpha*-lambda1.eqlin) + W(4,2)*(x2final-alpha*-lambda2.eqlin) ...
                   + W(4,3)*(x3final-alpha*-lambda3.eqlin) + W(4,4)*(x4final-alpha*-lambda4.eqlin);
        for i = 1:length(xf4(:,k+1))
            xf4(i,k+1) = max(xf4(i,k+1), xf4_min(i));
            xf4(i,k+1) = min(xf4(i,k+1), xf4_max(i));
        end
        
    end
    errors_total = [errors_total, errors];
end


%% ADMM for Consensus Optimization

iter = 700;
rhos = 0.1;
%rhos = [0.5, 0.1, 0.05, 0.01, 0.005];
errors_total = [];

for a = 1:length(rhos)
    
    % initialize lambdas
    y1 = zeros(4,iter);
    y2 = zeros(4,iter);
    y3 = zeros(4,iter);
    y4 = zeros(4,iter);
    subs1 = zeros(4,iter);
    subs2 = zeros(4,iter);
    subs3 = zeros(4,iter);
    subs4 = zeros(4,iter);
    errors = zeros(4,iter);
    rho = rhos(a);
    x_bar = zeros(4,1);

    for k = 1:iter
        % node 1 determines u1 that minimizes the objective and returns xf1
        cvx_clear
        cvx_begin quiet
            variable u1(2*N) 
            minimize ( (S1*u1+T1*x01)'*(S1*u1+T1*x01) + u1'*u1 + y1(:,k)'*(S1(end-3:end,:)*u1+T1(end-3:end,:)*x01-x_bar) + rho/2*(S1(end-3:end,:)*u1+T1(end-3:end,:)*x01-x_bar)'*(S1(end-3:end,:)*u1+T1(end-3:end,:)*x01-x_bar) )
            subject to
                norm(u1,Inf) <= umax/Tfinal;
        cvx_end

        xf1 = S1(end-3:end,:)*u1+T1(end-3:end,:)*x01;

        % node 2 
        cvx_clear
        cvx_begin quiet
            variable u2(2*N) 
            minimize ( (S2*u2+T2*x02)'*(S2*u2+T2*x02) + u2'*u2 + y2(:,k)'*(S2(end-3:end,:)*u2+T2(end-3:end,:)*x02-x_bar) + rho/2*(S2(end-3:end,:)*u2+T2(end-3:end,:)*x02-x_bar)'*(S2(end-3:end,:)*u2+T2(end-3:end,:)*x02-x_bar) )
            subject to
                norm(u2,Inf) <= umax/Tfinal;
        cvx_end 

        xf2 = S2(end-3:end,:)*u2+T2(end-3:end,:)*x02 ;

        % node 3 
        cvx_clear
        cvx_begin quiet
            variable u3(2*N) 
            minimize ((S3*u3+T3*x03)'*(S3*u3+T3*x03) + u3'*u3 + y3(:,k)'*(S3(end-3:end,:)*u3+T3(end-3:end,:)*x03-x_bar) + rho/2*(S3(end-3:end,:)*u3+T3(end-3:end,:)*x03-x_bar)'*(S3(end-3:end,:)*u3+T3(end-3:end,:)*x03-x_bar) )
            subject to
                norm(u3,Inf) <= umax/Tfinal;
        cvx_end 

        xf3 = S3(end-3:end,:)*u3+T3(end-3:end,:)*x03;

        % node 4
        cvx_clear
        cvx_begin quiet
            variable u4(2*N) 
            minimize ((S4*u4+T4*x04)'*(S4*u4+T4*x04) + u4'*u4 + y4(:,k)'*(S4(end-3:end,:)*u4+T4(end-3:end,:)*x04-x_bar) + rho/2*(S4(end-3:end,:)*u4+T4(end-3:end,:)*x04-x_bar)'*(S4(end-3:end,:)*u4+T4(end-3:end,:)*x04-x_bar) )
            subject to
                norm(u4,Inf) <= umax/Tfinal;
        cvx_end 

        xf4 = S4(end-3:end,:)*u4+T4(end-3:end,:)*x04 ;

        % store value of the errors
        subs1(:,k) = x_bar-xf1;
        subs2(:,k) = x_bar-xf2;
        subs3(:,k) = x_bar-xf3;
        subs4(:,k) = x_bar-xf4;
        errors(:,k) = [norm(subs1(:,k),1); norm(subs2(:,k),1); norm(subs3(:,k),1); norm(subs4(:,k),1)];
        
        % update x_bar (the average of all the final states)
        x_bar = 1/4*(xf1 + xf2 + xf3 + xf4);
        
        % update y variables
        y1(:,k+1) = y1(:,k) + rho*(xf1-x_bar);
        y2(:,k+1) = y2(:,k) + rho*(xf2-x_bar);
        y3(:,k+1) = y3(:,k) + rho*(xf3-x_bar);
        y4(:,k+1) = y4(:,k) + rho*(xf4-x_bar);
    end
    
    errors_total = [errors_total; errors];
end


%% Plot error sequences and the plane trajectories after convergence

figure(4)
plot((errors)') 
title('Error sequence per constraint')
legend('x_{final,1} - x_{final,2}', 'x_{final,2} - x_{final,3}', 'x_{final,3} - x_{final,4}', 'x_{final,4} - x_{final,1}')
% 'xf1 - x_bar', 'xf2 - x_bar', 'xf3 - x_bar', 'xf4 - x_bar'
xlabel('Iteration')
ylabel('Error')

x1 = S1*u1+T1*x01;
x2 = S2*u2+T2*x02;
x3 = S3*u3+T3*x03;
x4 = S4*u4+T4*x04;

dim_x = 4;
plt_x_dist = [];
plt_y_dist = [];

for i = 0:N-1
    plt_x_dist = [plt_x_dist; x1(i*dim_x+1),x2(i*dim_x+1),x3(i*dim_x+1),x4(i*dim_x+1)];
end
for i = 0:N-1
    plt_y_dist = [plt_y_dist; x1(i*dim_x+2),x2(i*dim_x+2),x3(i*dim_x+2),x4(i*dim_x+2)];
end

figure(5)
plot(plt_x_dist)
legend('plane 1', 'plane 2', 'plane 3', 'plane 4')
title('Position of the planes in the x-direction')
xlabel('Time step')
ylabel('Position')

figure(6)
plot(plt_y_dist)
legend('plane 1', 'plane 2', 'plane 3', 'plane 4')
title('Position of the planes in the y-direction')
xlabel('Time step')
ylabel('Position')

figure(7)
plot(plt_x_dist,plt_y_dist,'o-')
legend('plane 1', 'plane 2', 'plane 3', 'plane 4')
title('Position of the planes over time')
xlabel('Position x-direction')
ylabel('Position y-direction')


%% Plot error sequences for different step sizes/ step size update sequences

% plot of the error sequences per plane
figure(6)
plot(errors_total(4,:)')
hold on
plot(errors_total(8,:)')
plot(errors_total(12,:)')
plot(errors_total(16,:)')
plot(errors_total(20,:)')
title('Error sequences of plane 4 for different \rhos')
xlabel('Error')
ylabel('Iteration')
legend('rho = 0.5','rho = 0.1' ,'rho = 0.05','rho = 0.01','rho = 0.005')
% alternative legend: '\alpha(t) = 1/(10+t)','\alpha(t) = 0.1/sqrt(t)'

