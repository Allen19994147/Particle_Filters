%%% Particle filter implementation
clear all
clc
%%%%-----------Prerequisite-----------%%%%
x = [0;0];  %x as posi and velo
cycle = 50;  % repeatation
M = 100;    % #particles
uq = 0;     %mean of the measurement
q = 10;     %variance for measurement
uw = 0;     %mean of disturbance
w = 1;      %Variance of the disturbance
dt = 1;


%%%-----------------
z = [0;0];  %measurement
C = eye(2);
V = [1;1];
%z=[0];
%C = [1 0];
%V = 1;
%%%-----------------
A = [1 dt;0 1];
W = [0.5*dt*dt;dt];
R = w.*W*W';
Q = q.*V*V';

%%%--------Estimator--------------------%%%
x_est = [0;0];
z_est = [];

Var = 0.5;  % Initail state variance is zero

for i = 1:M
    for j = 1:2
        x_P(j,i) = x(j) + sqrt(Var) * randn; 
    end
end
x_P;    %%% Check the distribution of initial states
x_P_update = [0;0];
x_P_all = x_P;
%%%
for i=1:cycle
    i;
    x = [x,A*x(:,i)+normrnd(uw,w).*W(1)];  % Real states of robot
    z = [z,C*x(:,i+1)+normrnd(uq,q).*V];  % Measure states of robot
    
    z(:,i+1);       %it's latest measurement
    for j = 1:M
        x_P_update(1,j) = A(1,:)*x_P(:,j)+normrnd(uw,w).*W(1); % particle posi simulation
        x_P_update(2,j) = A(2,:)*x_P(:,j)+normrnd(uw,w).*W(2); % particle velo simulation
    
        %%%simu measurements should not include meazurement noises
        z_update(1,j) = C(1,:)*x_P_update(:,j);
        z_update(2,j) = C(2,:)*x_P_update(:,j);
        
        
        p_posi =  1/sqrt(2*pi*q) * exp(-(z(1,i+1) - z_update(1,j))^2/(2*q));
        p_velo =  1/sqrt(2*pi*q) * exp(-(z(2,i+1) - z_update(2,j))^2/(2*q));
        p_w(j) = sqrt(p_posi*p_velo);   % approx posibility of posi&velo
        
    end
    p_w = p_w./sum(p_w);       % weight normalization
    
    s = RandStream('mlfg6331_64'); 
    x_P(1,:) = randsample(x_P_update(1,:),M,true,p_w);
    x_P(2,:) = randsample(x_P_update(2,:),M,true,p_w);
    
    x_P_all= [x_P_all ; x_P];   % Track all the iter of particles
    
    x_est(1,i+1) = mean(x_P(1,:));
    x_est(2,i+1) = mean(x_P(2,:));
end

%%% plot
x;      %%%The true states
x_est;  %%%The estimated states(final results)
size(x_P_all);
% figure(1)
% for k = 1:cycle+1
%     subplot(1,cycle+1,k)
%     scatter(x_P_all(2*k-1,:),x_P_all(2*k,:),'b')
%     hold on
%     scatter(x(1,k),x(2,k),'r')
%     scatter(x(1,:),x(2,:),'r')
%     xlim([-15 15])
%     ylim([-8 8])
%     xlabel('position'); ylabel('velocity');
% end
%%%----------------------------
thesize = size(x_est);
i = 1:thesize(2);
figure(2);
clf
plot(i, x(1,:), '-.o', i, x_est(1,:), '-*','linewidth',3);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Iteration'); ylabel('position');
legend('True position', 'Particle estimate');

figure(3);
clf
plot(i, x(2,:),'-.o', i, x_est(2,:),'-*','linewidth',3);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Iteration'); ylabel('velocity');
legend('True velocity', 'Particle estimate');

%%% Square error estimation
x_err = [0;0];
for i = 1:cycle+1       %%%Calculate cumulative square errors
    x_err(1,1) = x_err(1,1) + (x(1,i)-x_est(1,i))^2;
    x_err(2,1) = x_err(2,1) + (x(2,i)-x_est(2,i))^2;
end
x_err = x_err./(cycle+1)


