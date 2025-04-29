function [X,P] = prediction_UKF(X,P,G,dt,alpha,beta,kappa,sigma_v,sigma_omega,sigma_ext)
    % sigma point parameters
    nx = length(X);
    lambda = alpha^2*(nx + kappa) - nx;
    
    % calculate weights for sigma points
    weights_m = zeros(1,2*nx+1); weights_c = weights_m;
    weights_m(1) = lambda/(nx + lambda);
    weights_c(1) = lambda/(nx + lambda) + (1 - alpha^2 + beta);
    weights_m(2:2*nx + 1) = 1/(2*(nx + lambda));
    weights_c(2:2*nx + 1) = 1/(2*(nx + lambda));

    %Calculate Sigma Points
    A = sqrt(nx + lambda) * chol(P)';
    sig_points = [zeros(size(X)) -A A];
    sig_points = sig_points + repmat(X, 1, size(sig_points, 2));
    
    % predict samples according to transition model
    for i = 1:2*nx+1
        sig_points(:,i) = motion_model(sig_points(:,i),dt);
    end

    % calculate process noise covariance matrix
    G(1,1) = G(1,1)*cos(X(5));
    G(2,1) = G(2,1)*sin(X(5));
    Q = G*diag([sigma_v sigma_omega sigma_ext*ones(1,nx-6)].^2)*G';
    
    % calculate predicted mean
    X = sum(weights_m.*sig_points,2);
    P = (weights_c.*(sig_points - X))*(sig_points - X)';
    P = P + Q;
end