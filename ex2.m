%%%HW3 Kalman filter (FULL) implementation
clear all
clc
clf
%%%%-----------Prerequisite-----------%%%%
x = [0;0];  %x as posi and vel
w = [0;0];  %model accleration
v = [0;0];  %measure noise
dt = 1;
A = [1 dt;0 1];
W = [0.5*dt*dt;dt];%Disturbance matrix
%%%-----------------
% z = [0;0];  %measurement
% C = eye(2);
% V = [1;1];
z=[0];
C = [1 0];
V = 1;
%%%-----------------   
wu = 0;     %mean of the acc/input
w = 1;      %variance of disturbance
uq = 0;     %mean of the measurement
q = 10;     %variance for measurement
R = w.*W*W';
Q = q.*V*V';

cycle = 50;
%%%--------Estimator--------------------%%%
x_est = [0;0];
z_est = [];

K_all = [];
%%%
x_cov = [0 0;0 0];  %Initial covariance is zero
Cov = x_cov;    %Track all the valus over time
Cov_bar = x_cov;
for i=1:cycle
    i
    theX = A*x(:,i)+normrnd(wu,w).*W;   %Real state of robots
    x = [x,theX];
    z = [z,C*x(:,i+1)+normrnd(uq,q).*V];  % Measure states of robot
    
    %%%%-------------------
    %x_est_bar = A*x_est(:,i);   %state estimation without noise
    x_est_bar = A*x_est(:,i) + normrnd(wu,w).*W;    %state estimation with noise
%     if(i==5)
%         z(:,6) = 5;
%         x(:,6) = [0;0];
%     end
    x_cov_bar = A*x_cov*A' + R
    Cov_bar = [Cov_bar, x_cov_bar];
    K = x_cov_bar*C'*inv(Q + C*x_cov_bar*C');
    K_all = [K_all, K];
    
    %Estimated states with correction
    x_est = [x_est, x_est_bar+K*(z(:,i+1)-C*x_est_bar)];      
    x_cov = (eye(2) - K*C)*x_cov_bar      
    Cov = [Cov, x_cov];
end
% K_all
% x
% z
% x_est
Cov
% Cov_bar

%%%%%%%%%%plot real and est here%%%%%%%%%%
thesize = size(x_est);
i = 1:thesize(2);

figure(1);
clf
plot(i, x(1,:), '-.o', i, x_est(1,:), '-*','linewidth',3);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('Iteration'); ylabel('position');
legend('True position', 'KF estimate');

%%% Square error estimation
x_err = [0;0];
for i = 1:cycle+1       %%%Calculate cumulative square errors
    x_err(1,1) = x_err(1,1) + (x(1,i)-x_est(1,i))^2;
    %x_err(2,1) = x_err(2,1) + (x(2,i)-x_est(2,i))^2;
end
x_err = x_err./(cycle+1)


% %%%----------------error eclipse--------------------------------
figure(2)
for i=1:cycle
%繪製error ellipse
% Calculate the eigenvectors and eigenvalues
% data = 
i;
covariance = [Cov(:,2*i+1) Cov(:,2*i+2)];
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the coordinates of the data mean
% avg = mean(data);

% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi);
phi = angle;
X0= x_est(1,i+1);
Y0= x_est(2,i+1);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
if i==5
    
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-')
end
hold on;

% Set the axis labels
hXLabel = xlabel('x');
hYLabel = ylabel('y');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:cycle
%繪製error ellipse
% Calculate the eigenvectors and eigenvalues
% data = 
i;
covariance = [Cov_bar(:,2*i+1) Cov_bar(:,2*i+2)];
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the coordinates of the data mean
% avg = mean(data);

% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi);
phi = angle;
X0= x_est(1,i+1);
Y0= x_est(2,i+1);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
if i==5
    
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'*')
end
hold on;

% Set the axis labels
hXLabel = xlabel('x');
hYLabel = ylabel('y');
legend("after measurement","before measurement")
end
hold off

