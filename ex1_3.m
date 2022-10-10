%%%HW3 Kalman filter (pre) implementation
clear all
clc
clf
%%%%-----------Prerequisite-----------%%%%
x = [0;0];  %x as posi and velo
u = [0;0];  %u as acceleration
w = [0;0];  %model disturbance
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
w = 1;      %variance for acc/input
uq = 0;     %mean of the measurement
q = 10;     %variance for measurement
R = w.*W*W';
Q = q.*V*V';
cycle = 5;
%%%--------Estimator--------------------%%%
x_est = [0;0];
z_est = [];
K_all = [];
%%%
x_cov = [0 0;0 0];  %Initial covariance is zero
Cov = x_cov;    %Track all the valus over time
for i=1:cycle
    i
    theX = A*x(:,i)+normrnd(wu,w).*W;   % Real states of robot
    x = [x,theX];   %Record all of them
    z = [z,C*x(:,i+1)+normrnd(uq,q).*V];  % Measure states of robot
    %%%%-------------------
    x_est_bar = A*x_est(:,i);
    x_est = [x_est,x_est_bar];   % states estimation
    x_cov = A*x_cov*A' + R; % current covariance
    Cov = [Cov,x_cov];  % Record all variance
end
%%% followings are only for demonstration purpose
Cov(:,1:2) = [];
Cov      %%%   Remove semicolon to see Cov
Cov = [[0 0;0 0], Cov];
x
x_est


%%%----------------error eclipse--------------------------------
MIN_eiganVal = [];
MAX_eiganVal = [];
for i=1:cycle
%繪製error ellipse
% Calculate the eigenvectors and eigenvalues
i;
covariance = [Cov(:,2*i+1) Cov(:,2*i+2)]
[eigenvec, eigenval ] = eig(covariance);

% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));
MAX_eiganVal = [MAX_eiganVal,largest_eigenval];
% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end
MIN_eiganVal = [MIN_eiganVal,smallest_eigenval];
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
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'-')
hold on;

% Set the axis labels
hXLabel = xlabel('position');
hYLabel = ylabel('velocity');
end
MIN_eiganVal
MAX_eiganVal