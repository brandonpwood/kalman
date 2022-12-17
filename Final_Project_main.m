close all
clear
clc

% ECE 5250/4250 Final project

% Specify the Boundary box specification
B_x = 20;                   % for open boundary case, set this value to inf
                            % for closed boundary case, specify the max x value, 
B_y = 20;                   % for open boundary case, set this value to inf
                            % for closed boundary case, specify the max y value,

% Specify the model parameters, A, C

% s[n+1] = A*s[n] + x[n] + u[n]    u[n]~N(0,R[n])
% y[n] = C*s[n] + w[n]             w[n]~N(0,Q[n])
A = [1 0 1 0;
     0 1 0 1;
     0 0 1 0;
     0 0 0 1];
C = [1 0 0 0;
     0 1 0 0];

% Specify the Number of targets
num_tar = 3;

% Specify the process noise covariance matrix, R for the individual targets
% if the number of targets > 3, specify additional diagonal matrices R(:,:,i)
R(:,:,1) = diag([0.001, 0.001, 0.001, 0.001]);  
R(:,:,2) = diag([0.001, 0.001, 0.001, 0.001]);
R(:,:,3) = diag([0.001, 0.001, 0.001, 0.001]);

% Specify the measurement noise covariance matrix, Qfor the individual targets
% if the number of targets > 3, specify additional diagonal matrices Q(:,:,i)
Q(:,:,1) = diag([0.005, 0.005]);  
Q(:,:,2) = diag([0.05, 0.05]);
Q(:,:,3) = diag([0.5, 0.5]);

% Specify the initial state of the targets in the following order 
% x position, y position, x velocity, y velcoity
% if the number of targets > 3, specify additional intial conditions
s0 = [ 0  5  0.1  0.05;            % Target 1, x position, y position, x velocity, y velcoity
       0  5   0.1  0.05;               % Target 2, x position, y position, x velocity, y velcoity
       10  0  -0.1  0.05];             % Target 3, x position, y position, x velocity, y velcoity

% Specify the duration of simulation, and and the time interval, dt
duration = 250; 

% n is the duration. You don't have to do anything here
n = duration;              % number of time steps

%% Data generation

% input : [B_x, B_y, A, C, R, Q, num_tar, s0, n]
          % Bx, By are the x and y boundaries.
          % A: 4 X 4 state transition matrix.
          % C: 2 X 4 observation matrix.
          % R: 4 X 4 x number of targets, 3D process noise covariance matrix.
          % Q: 2 X 24 x number of targets, observation noise covariance matrix.
          % num_tar: number of targets.
          % s0: tar X 4 initial states.
          % n: duration.

% output: [store_st_mm, x]
          % store_st_mm: 3D matrix (dimesnion: 6 x number of time steps x number of targets.)
            % Data generated for different targets is stored in different
            % layers of this matrix. The individual layer consists of six rows. 
            % The first four rows are the states, while the next two rows are 
            % the noisy measurements. 
          % x: 3D Matrix (dimension: 4 x number of time steps x number of targets)
            % represents the control inputs for different targets.

[Store_st_mm, x] = data_generation(B_x, B_y, A, C, R, Q, num_tar, s0, n);


%% Kalman Filter for a single target

% You should develope a single Kalman filter function for all the targets with the 
% following as the inputs and outputs of your Kalman filter function

% Kalman filter function
% input : [A, R, C, Q, Store_st_mm, x, num_tar, n]
          % A: 4 X 4 state transition matrix.
          % C: 2 X 4 observation matrix.
          % R: 4 X 4 x number of targets, 3D process noise covariance matrix.
          % Q: 2 X 24 x number of targets, observation noise covariance matrix.
          % store_st_mm: 3D matrix (dimesnion: 6 x number of time steps x number of targets.)
            % Data generated for different targets is stored in different
            % layers of this matrix. The individual layer consists of six rows. 
            % The first four rows are the states, while the next two rows are 
            % the noisy measurements.
          % x: 3D Matrix (dimension: 4 x number of time steps x number of targets)
            % represents the control inputs for different targets.
          % num_tar: number of targets.
          % n: duration

% output: [s_hat and sig_hat, K]
          % s_hat: 2D matrix (dimesnion: (4 num_tar) x number of time steps)
            % The columns of s_hat should contain the state estimates of the targets. 
            % The first four states in a column should correspond to the state estimate of 
            % the first target. The next four states should correspond to the state estimates
            % for the second target and so on..
          % sig_hat: 3D matrix (dimesnion: (4 num_tar) x (4 num_tar) x number of time steps) 
            % The 2D layers of this 3D matrix should represent the covariance matrices 
            % for individual time steps. Block diagonal elements at a particular time 
            % step, i.e., in a particular layer, should correspond to the covariance 
            % matrices of different targets.
          % K: 3D matrix (dimension: (4 num_tar) x (2 num_tar) x number of time steps)
            % The 2D layers of this matrix should represent the Kalman Gain matrix for 
            % the individual time steps. The blocks within a layer should should correspond
            % to the Kalman Gain matrices of different targets.
 
% specify the target number for Kalman filtering 
[s_hat, sig_hat, K] = Kalman_filt(A, R, C, Q, Store_st_mm, x, num_tar, n);


%% Monte Carlo part

% input : [A, R, C, Q, Store_st_mm, x, num_tar, n, num_sim]
          % A: 4 X 4 state transition matrix.
          % C: 2 X 4 observation matrix.
          % R: 4 X 4 x number of targets, 3D process noise covariance matrix.
          % Q: 2 X 24 x number of targets, observation noise covariance matrix.
          % store_st_mm: 3D matrix (dimesnion: 6 x number of time steps x number of targets.)
            % Data generated for different targets is stored in different
            % layers of this matrix. The individual layer consists of six rows. 
            % The first four rows are the states, while the next two rows are 
            % the noisy measurements.
          % x: 3D Matrix (dimension: 4 x number of time steps x number of targets)
            % represents the control inputs for different targets.
          % num_tar: number of targets.
          % n: duration.
          % num_sim: number of samples to be used in the Monte Carlo simulation.

% Output: [MMSE_Monte_Carlo]
          % MMSE_Monte_Carlo: 2D matrix (dimension: (num_tar) x number of time steps)
            % MMSEs of the states of the targets at a specific time step are given by 
            % the columns of this matrix. The rows of this amtrix correspond to MMSEs from 
            % Monte Carlo method for different targets.
        
num_sim = 500;

MMSE_Monte_Carlo = monte_carlo(A, R, C, Q, Store_st_mm, x, num_tar, n, num_sim);

%% Function for plotting p1.1, p1.2, p1.4, p1.5

% Input: [B_x, B_y, num_tar, Store_st_mm, s_hat, sig_hat, K, MMSE_Monte_Carlo, n, t]
          % Bx, By are the x and y boundaries.
          % num_tar: number of targets.
          % Store_st_mm: Output from data generation
          % s_hat, sig_hat, K: Outputs from Kalman Filter
          % MMSE_Monte_Carlo: output from Monte Carlo Simulations
          % n: duration.
          % t: time vector,[1,2,...250]

% Output: [Plots corresponding to p1.1, p1.2, p1.4, p1.5]

t = [1:duration];       
plotting(B_x, B_y, num_tar, Store_st_mm, s_hat, sig_hat, K, MMSE_Monte_Carlo, n, t)
