% function for data generation

% input : [B_x, B_y, A, C, R, Q, num_tar, s0, n]
          % Bx, By are the x and y boundaries.
          % A: 4 X 4 state transition matrix.
          % C: 2 X 4 observation matrix.
          % R: 4 X 4 x number of targets, 3D process noise covariance
          % matrix.
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

function [Store_st_mm, x] = data_generation(B_x, B_y, A, C, R, Q, num_tar, s0, n)

    for tar = 1:num_tar 
        % process noise, 0 mean Gaussian
        u(:,:,tar) = [randn(n,1)*sqrt(R(1, 1, tar))... 
                      randn(n,1)*sqrt(R(2, 2, tar))...
                      randn(n,1)*sqrt(R(3, 3, tar))...
                      randn(n,1)*sqrt(R(4, 4, tar))];    

        % measurement noise 0 mean Gaussian
        w(:,:,tar) = [randn(n,1)*sqrt(Q(1,1, tar))... 
                      randn(n,1)*sqrt(Q(2,2, tar))];  

    end

    % Simulating all the targets from 1 to num_tar
    Store_st_mm = zeros(2, n, num_tar);

    % control input
    x_c = zeros(4,n);  
    x = zeros(4,n,num_tar);

    for tar = 1:num_tar
    
        s = s0(tar,:)';
        y = zeros(2, n);

        for k=1:(n-1),
        
            s(:,k+1) = A*s(:,k) + x_c(:,k) + [u(k,:,tar)]';         % with process noise
            y(:,k+1) = C*s(:,k+1) + [w(k+1,:, tar)]';
    
            %If the target is outside the boundary
            if s(1,k+1) > B_x                      
                x_c(1,k) = -2 * ((s(1,k)+s(3,k))-B_x);
                x_c(3,k) = -2 * s(3,k);
                s(:,k+1) = A * s(:,k) + x_c(:,k) + [u(k,:,tar)]';      % with process noise
                y(:,k+1) = C * s(:,k+1) + [w(k+1,:, tar)]';
            end
        
            if s(2,k+1) > B_y  
                x_c(2,k) = -2 * ((s(2,k)+s(4,k))-B_y);
                x_c(4,k) = -2 * s(4,k);
                s(:,k+1) = A*s(:,k) + x_c(:,k) + [u(k,:, tar)]';      % with process noise
                y(:,k+1) = C*s(:,k+1) + [w(k+1,:, tar)]';
            end
    
            if s(1,k+1) < 0
                x_c(1,k) = -2 * ((s(1,k)+s(3,k)));
                x_c(3,k) = -2 * s(3,k);
                s(:,k+1) = A*s(:,k) + x_c(:,k) + [u(k,:, tar)]';      % with process noise
                y(:,k+1) = C*s(:,k+1) + [w(k+1,:, tar)]';
            end

            if s(2,k+1) < 0
                x_c(2,k) = -2 * (s(2,k)+s(4,k));
                x_c(4,k) = -2 * s(4,k);
                s(:,k+1) = A*s(:,k) + x_c(:,k) + [u(k,:, tar)]';      % with process noise
                y(:,k+1) = C*s(:,k+1) + [w(k+1,:, tar)]';
            end
    
        end

    % Storing the xtrue and z values from different targets
    Store_st_mm(1:4,:, tar) = s;        % First four rows are states
    Store_st_mm(5:6,:, tar) = y;        % Next two rows are the measurements

    x(:,:,tar) = x_c;

    end
    
end