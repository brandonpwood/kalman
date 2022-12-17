% Monte Carlo Simulation function

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
          % Store_st_mm: see above for details.
          % num_sim: number of samples to be used in the Monte Carlo simulation.

% Output: [MMSE_Monte_Carlo]
          % MMSE_Monte_Carlo: 2D matrix (dimension: (num_tar) x number of time steps)
            % MMSEs of the states of the targets at a specific time step are given by 
            % the columns of this matrix. The rows of this amtrix correspond to MMSEs from 
            % Monte Carlo method for different targets.

function MMSE_Monte_Carlo = monte_carlo(A, R, C, Q, Store_st_mm, x, num_tar, n, num_sim)

    %sq_error_mat = zeros(4,n,num_sim);
    Store_st_mm_c = Store_st_mm(1:4,:,1);
    for tar = 2 : num_tar
        Store_st_mm_c(end+1:end+4, :) = Store_st_mm(1:4,:,tar);
    end

    for j = 1:(num_sim)

        [shatm, sighatm, K] = Kalman_filt(A, R, C, Q, Store_st_mm, x, num_tar, n);

        % Mean square error from Monte Carlo Method
        sq_error_mat(:,:,j) = (Store_st_mm_c - shatm).^2;

    end

;

    all_error = mean(sq_error_mat,3);
    for tar = 1: num_tar
        MMSE_Monte_Carlo(tar,:) = sum(all_error((4*tar)-3:4*tar,:));
    end

end