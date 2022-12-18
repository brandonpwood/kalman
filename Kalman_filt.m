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


function [s_hat, sig_hat, K] = Kalman_filt(A, R, C, Q, store_st_mm, x, num_tar, n)
    % Initialize return functions to zero
    s_hat = zeros(4*num_tar, n);
    sig_hat = zeros( 4*num_tar, 4*num_tar, n);
    K = zeros(4*num_tar, 2*num_tar, n);
    
    % Keep track of predictions
    s_preds = zeros(4, num_tar, n); 
    y_preds = zeros(2, num_tar, n);
    sig_preds = zeros(4, 4, num_tar, n);
    for i = [1:n]
        if i == 1
            %Intialize slices for this point in time 
            K_acc = zeros(4*num_tar, 2*num_tar); 
            s_acc = [];
            sigma_acc = zeros(4*num_tar, 4*num_tar);
            
            % Iterate over targets        
            for tar = [1:num_tar]
                % Initial estimation
                s_0 = [1.25;2.25;.25;3.25];   
                y_0 = C*s_0;
                sig_0 = 10*ones(4, 4);
                
                % Recursive updates
                q = Q(:, :, tar);
                r = R(:, :, tar);
                sigma_j = sig_0;
                
                % Update K
                k_i = sigma_j*(C.')*((C*sigma_j *(C.') + q)^-1);
                s = size(k_i);
                K_acc = insert_diagonal(K_acc, (tar-1)*s(1), (tar-1)*s(2), k_i); % append to acculumator
                
                % Estimate State
                s_prev = s_0;
                y_prev = y_0;
                y_meas = store_st_mm(5:6, i, tar); % Extract measurements from input data
                s_est = s_prev + k_i*(y_meas - y_prev); % new best estimate of s
                s_acc = [s_acc; s_est];
                
                % Compute covariance
                sigma_est = sigma_j - k_i*C*sigma_j;
                s_sigma = size(sigma_est);
                sigma_acc = insert_diagonal (sigma_acc, (tar-1)*s_sigma(1), (tar-1)*s_sigma(2), sigma_est);
                
                
                % Compute new predictions
                s_pred = A*s_est;
                y_pred = C*s_pred;               
                sigma_pred = A*sigma_est*(A.') + r;
                
                s_preds (:, tar, i) = s_pred;
                y_preds(:, tar,i) = y_pred;
                sig_preds(:, :, tar, i) = sigma_pred;
            end
            % Reformat Arrays Across Targets            
            K(:, :, i) = K_acc; 
            s_hat (:, i) = s_acc;
            sig_hat(:, :, i) = sigma_acc;
   
        else
            %sigma = sig_hat(:, :, i-1);
            % Make arrays to store results for each target (state estimates)
            K_acc = zeros(4*num_tar, 2*num_tar); 
            s_acc = [];
            
            sigma_acc = zeros(4*num_tar, 4*num_tar);
            
            
            % Iterate over targets
            for tar = [1:num_tar]
                % Recursive updates
                q = Q(:, :, tar);
                r = R(:, :, tar);
                %slice_x = 1 + (tar-1)*4;
                %slice_y = slice_x + 3;
                %sigma_j = sigma(slice_x : slice_y, slice_x : slice_y);
                sigma_j = sig_preds(:, :, tar, i-1);
                % Update K
                k_i = sigma_j*(C.')* (inv (C*sigma_j *(C.') + q) );

                s = size(k_i);
                K_acc = insert_diagonal(K_acc, (tar-1)*s(1), (tar-1)*s(2), k_i); % append to acculumator
                
                % Estimate State
                s_prev = s_preds(:, tar, i-1);
                y_prev = y_preds(:, tar, i-1);
                y_meas = store_st_mm(5:6, i, tar); % Extract measurements from input data
                s_est = s_prev + k_i*(y_meas - y_prev); % new best estimate of s
                s_acc = [s_acc; s_est];
                
                % Compute covariance
                sigma_est = sigma_j - k_i*C*sigma_j;
                s_sigma = size(sigma_est);
                sigma_acc = insert_diagonal (sigma_acc, (tar-1)*s_sigma(1), (tar-1)*s_sigma(2), sigma_est);
                
                
                % Compute new predictions
                s_pred = A*s_est;
                y_pred = C*s_pred;               
                sigma_pred = A*sigma_est*(A.') + r;
                
                s_preds (:, tar, i) = s_pred;
                y_preds(:, tar,i) = y_pred;
                sig_preds(:, :, tar, i) = sigma_pred;
            end
            % Reformat Arrays Across Targets            
            K(:, :, i) = K_acc;
            K(:, :, i);
            s_hat (:, i) = s_acc;
            sig_hat(:, :, i) = sigma_acc;
            
        end
    end
    
end   