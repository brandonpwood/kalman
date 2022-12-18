close all
clear
clc


B_x = 20;                   % for open boundary case, set this value to inf
                            % for closed boundary case, specify the max x value, 
B_y = 20;                   % for open boundary case, set this value to inf
                            % for closed boundary case, specify the max y value,

A = [1 0 1 0;
     0 1 0 1;
     0 0 1 0;
     0 0 0 1];
C = [1 0 0 0;
     0 1 0 0];

% Make a few exrta targets to tets different covariances at once
num_tar = 3;

R(:,:,1) = diag([0.001, 0.001, 0.001, 0.001]);  
R(:,:,2) = diag([0.001, 0.001, 0.001, 0.001]);
R(:,:,3) = diag([0.001, 0.001, 0.001, 0.001]);

Q(:,:,1) = diag([0.05, 0.05]);  
Q(:,:,2) = diag([0.05, 0.05]);
Q(:,:,3) = diag([0.05, 0.05]);


s0 = [ 0  5  0.1  0.05;
10 0  -0.1  0.05;
15  15  -0.1  -0.05]; 

duration = 250; 
n = duration; 

[Store_st_mm, x] = data_generation(B_x, B_y, A, C, R, Q, num_tar, s0, n);

[s_hat, sig_hat, K] = Kalman_filt(A, R, C, Q, Store_st_mm, x, num_tar, n);

num_sim = 500;
MMSE_Monte_Carlo = monte_carlo(A, R, C, Q, Store_st_mm, x, num_tar, n, num_sim);


t = [1:duration];
plotting(B_x, B_y, num_tar, Store_st_mm, s_hat, sig_hat, K, MMSE_Monte_Carlo, n, t)