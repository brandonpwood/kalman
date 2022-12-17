% Function for plotting p1.1, p1.2, p1.4, p1.5

% Input: [B_x, B_y, num_tar, Store_st_mm, s_hat, sig_hat, K, MMSE_Monte_Carlo, n, t]
          % Bx, By are the x and y boundaries.
          % num_tar: number of targets.
          % Store_st_mm: Output from data generation
          % s_hat, sig_hat, K: Outputs from Kalman Filter
          % MMSE_Monte_Carlo: output from Monte Carlo Simulations
          % n: duration.
          % t: time vector,[1,2,...250]

% Output: [Plots corresponding to p1.1, p1.2, p1.4, p1.5]

function plotting(B_x, B_y, num_tar, Store_st_mm, s_hat, sig_hat, K, MMSE_Monte_Carlo, n, t)

    figure(1)
    set(gcf, 'position', [800 100 1000 1000])
    if num_tar == 1 && B_x == inf
        sgtitle('Single Target trajectory, No boundary')
    elseif num_tar == 1 && B_x ~= inf
        sgtitle('Single Target trajectory, With boundary')
    elseif num_tar ~= 1 && B_x == inf
        sgtitle('Multi Target trajectories, No boundary')
    else num_tar ~= 1 && B_x ~= inf
        sgtitle('Multi Target trajectories, With boundary')
    end

    hold on
    col = [0,1,0;
           1,0,0;
           0,0,1];
    for k = 1:num_tar
        plot(Store_st_mm(1,:,k), Store_st_mm(2,:,k),'linewidth',2, 'Color', col(k,:));
    end
    plot([0, B_x, B_x, 0, 0],[0, 0, B_y, B_y, 0],'k-', 'linewidth', 2)

    if num_tar == 1
        legend('Target 1 trajectory','Boundary','location','best');
    else
        legend('Target 1 trajectory','Target 2 trajectory', 'Target 3 trajectory',...
            'Boundary','location','best');
    end
    xlabel('x position, [m]');
    ylabel('y position, [m]');
    grid on

    for tar = 1:num_tar
    
        figure(tar+1)
        set(gcf, 'position', [800 100 1000 1000])
        sgtitle(sprintf('Target - %d simulation results',tar))

        subplot(2,2,1)
        K11 = reshape(K(((4*tar)-3),((2*tar)-1),:),[1, n]);
        K12 = reshape(K(((4*tar)-3),2*tar,:),[1, n]);
        K21 = reshape(K(((4*tar)-2),((2*tar)-1),:),[1, n]);
        K22 = reshape(K(((4*tar)-2),2*tar,:),[1, n]);
        plot(t, K11,'g-','linewidth',2)
        hold on
        plot(t, K12,'b-' ,'linewidth',2)
        plot(t, K21,'c--','linewidth',2)
        plot(t, K22,'r--' ,'linewidth',2)
        legend('K_{11}', 'K_{12}','K_{21}', 'K_{22}')
        title('P 1.1: Kalman Gain vs Time')
        xlabel('time, [s]')
        ylabel('Kalman Gain')
        ylim([0,1])
        grid on

        % MMSE vs Time, from Kalman Filter
        subplot(2,2,2)
        sx = reshape(sig_hat(((4*tar)-3),((4*tar)-3),:),[1, n]);
        sy = reshape(sig_hat(((4*tar)-2),((4*tar)-2),:),[1, n]);
        svx = reshape(sig_hat(((4*tar)-1),((4*tar)-1),:),[1, n]);
        svy = reshape(sig_hat(4*tar,4*tar,:),[1, n]);
        plot(t, sx,'r-','linewidth',2)
        hold on
        plot(t, sy,'g--' ,'linewidth',2)
        plot(t, svx,'b-','linewidth',2)
        plot(t, svy,'m--' ,'linewidth',2)
        legend('\Sigma_{11}', '\Sigma_{22}','\Sigma_{33}', '\Sigma_{44}')
        title('P 1.2: MMSE vs Time, from Kalman Filter')
        xlabel('time, [s]')
        ylabel('MMSE')
        grid on
    
        % Actual trajectory and estimated trajectory
        subplot(2,2,3)
        hold on
        plot(Store_st_mm(1,:,tar), Store_st_mm(2,:,tar),'-','linewidth',2, 'Color', col(tar,:))
        plot(s_hat((4*tar)-3,:), s_hat((4*tar)-2,:),'-.','linewidth',2, 'Color', col(tar,:));
        plot([0, B_x, B_x, 0, 0],[0, 0, B_y, B_y, 0],'k-', 'linewidth', 2)

        legend('Actual State Trajectory','Estimated State Trajectory','location','best')
        title('P 1.4: Actual trajectory and estimated trajectory')
        xlabel('x position, [m]')
        ylabel('y position, [m]')
        grid on
        
        s_k = sx + sy + svx +svy;
    
        subplot(2,2,4)
        hold on
        plot(t, s_k,'g-','linewidth',2);
        plot(t, MMSE_Monte_Carlo(tar,:),'g-.','linewidth',2);
        legend('\Sigma','\Sigma Monte carlo', 'location', 'best')
        title('P 1.5: \Sigma vs \Sigma Monte carlo')
        xlabel('time, [s]')
        ylabel('MMSE')
        grid on

    end

end