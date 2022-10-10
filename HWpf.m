%%%HW3 Particle filter implementation for HW
clear all
clc
%%%%-----------Prerequisite-----------%%%%
x = [0;0];  %x as posi and velo
cycle = 5;  % repeatation
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
Var = 0;  % Initail state variance is zero
for i = 1:M
    for j = 1:2
        x_P(j,i) = x(j) + sqrt(Var) * randn; % initial partical distribution
    end
end
x_P_update = [0;0];
x_P_all = x_P;  %Record all particles over time
%%%
for i=1:cycle
    i;
    x = [x,A*x(:,i)+normrnd(uw,w).*W(1)];  % Real states of robot
    z = [z,C*x(:,i+1)+normrnd(uq,q).*V];  % Measure states of robot
    z(:,i+1);       %it's latest measurement
    for j = 1:M
        x_P_update(1,j) = A(1,:)*x_P(:,j)+normrnd(uw,w).*W(1); % particle posi simulation
        x_P_update(2,j) = A(2,:)*x_P(:,j)+normrnd(uw,w).*W(2); % particle velo simulation
        z_update(1,j) = C(1,:)*x_P_update(:,j); % position measurement update
        z_update(2,j) = C(2,:)*x_P_update(:,j); % velocity measurement update
        p_posi =  1/sqrt(2*pi*q) * exp(-(z(1,i+1) - z_update(1,j))^2/(2*q));    % position probability 
        p_velo =  1/sqrt(2*pi*q) * exp(-(z(2,i+1) - z_update(2,j))^2/(2*q));    % velocity probability 
        p_wp(j) = p_posi;   % position weight
        p_wv(j) = p_velo;   % velocity weight
    end
    p_wp = p_wp./sum(p_wp);       % weight normalization
    p_wv = p_wv./sum(p_wv);       % weight normalization
    
    s = RandStream('mlfg6331_64'); 
    x_P(1,:) = randsample(x_P_update(1,:),M,true,p_wp); %position resampling
    x_P(2,:) = randsample(x_P_update(2,:),M,true,p_wv); %velocity resampling
    
    x_P_all= [x_P_all ; x_P];   % Track all the iter of particles
    
    x_est(1,i+1) = mean(x_P(1,:));  %position mean 
    x_est(2,i+1) = mean(x_P(2,:));  %velocity mean 
end

%%% plot
x
x_est
size(x_P_all);
figure(1)
for k = 1:cycle+1
    subplot(1,cycle+1,k)
    scatter(x_P_all(2*k-1,:),x_P_all(2*k,:),'b')
    hold on
    scatter(x(1,k),x(2,k),'r')
    xlim([-15 15])
    ylim([-8 8])
    xlabel('position'); ylabel('velocity');
end
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

