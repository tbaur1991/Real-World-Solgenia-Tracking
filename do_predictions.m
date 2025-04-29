if strcmp(methods,'own')
    % predict elliptic cylinder filter
    [X_cyl(:,k),P_cyl(:,:,k)] = prediction_UKF(X_cyl(:,k-1),P_cyl(:,:,k-1),G_cyl,dt,alpha_UKF,beta_UKF,kappa,sigma_v,sigma_omega,sigma_ext);
    % predict FCDS filter
    if strcmp(source,'radial')
        if strcmp(filter,'GAM')
            [X_fcds(:,k),P_fcds(:,:,k)] = prediction_UKF(X_fcds(:,k-1),P_fcds(:,:,k-1),G_fcds,dt,alpha_UKF,beta_UKF,kappa,sigma_v,sigma_omega,sigma_ext);
        elseif strcmp(filter,'ERHM')
            % build full FCDS system state distribution
            X = [X_fcds(1:2,k-1); X_line(1,k-1); X_fcds(3:5,k-1); X_line(2,k-1); X_fcds(6:end,k-1)];
            P = zeros(size(X_fcds,1)+2,size(X_fcds,1)+2);
            P(1:2,1:2) = P_fcds(1:2,1:2,k-1); P(1:2,4:6) = P_fcds(1:2,3:5,k-1); P(1:2,8:end) = P_fcds(1:2,6:end,k-1);
            P(3,3) = P_line(1,1,k-1); P(3,7) = P_line(1,2,k-1);
            P(4:6,1:2) = P_fcds(3:5,1:2,k-1); P(4:6,4:6) = P_fcds(3:5,3:5,k-1); P(4:6,8:end) = P_fcds(3:5,6:end,k-1);
            P(7,3) = P_line(2,1,k-1); P(7,7) = P_line(2,2,k-1);
            P(8:end,1:2) = P_fcds(6:end,1:2,k-1); P(8:end,4:6) = P_fcds(6:end,3:5,k-1); P(8:end,8:end) = P_fcds(6:end,6:end,k-1);
            % do prediction step
            [X,P] = prediction_UKF(X,P,G_fcds,dt,alpha_UKF,beta_UKF,kappa,sigma_v,sigma_omega,sigma_ext);
            % save predicted system state distribution parameters
            X_fcds(:,k) = [X(1:2); X(4:6); X(8:end)]; 
            X_line(:,k) = [X(3); X(7)];
            P_fcds(1:2,1:2,k) = P(1:2,1:2); P_fcds(1:2,3:5,k) = P(1:2,4:6); P_fcds(1:2,6:end,k) = P(1:2,8:end);
            P_fcds(3:5,1:2,k) = P(4:6,1:2); P_fcds(3:5,3:5,k) = P(4:6,4:6); P_fcds(3:5,6:end,k) = P(4:6,8:end);
            P_fcds(6:end,1:2,k) = P(8:end,1:2); P_fcds(6:end,3:5,k) = P(8:end,4:6); P_fcds(6:end,6:end,k) = P(8:end,8:end);
            P_line(:,:,k) = [P(3,3) P(3,7); P(7,3) P(7,7)];
        end
    end
    % predict superellipse filter
    [X_super(:,k),P_super(:,:,k)] = prediction_UKF(X_super(:,k-1),P_super(:,:,k-1),G_super,dt,alpha_UKF,beta_UKF,kappa,sigma_v,sigma_omega,sigma_ext);
elseif strcmp(methods,'comp')
    % predict SDFS filter
    [X_sdfs(:,k),P_sdfs(:,:,k)] = prediction_UKF(X_sdfs(:,k-1),P_sdfs(:,:,k-1),G_sdfs,dt,alpha_UKF,beta_UKF,kappa,sigma_v,sigma_omega,sigma_ext);
    % predict SSH filter
    [X_sh(:,k),P_sh(:,:,k)] = prediction_UKF(X_sh(:,k-1),P_sh(:,:,k-1),G_sh,dt,alpha_UKF,beta_UKF,kappa,sigma_v,sigma_omega,sigma_ext);
    % predict 3D GP filter
    [X_gp(:,k),quat_gp(:,k),P_gp(:,:,k)] = prediction_3DGP(X_gp(:,k-1),quat_gp(:,k-1),P_gp(:,:,k-1),F,Q,lambda_gp,dt,stdAngVel);
end