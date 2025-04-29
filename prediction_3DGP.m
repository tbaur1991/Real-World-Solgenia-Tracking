function [predState,predQuat,predStateCov] = prediction_3DGP(prevEstState,prevEstQuat,prevEstStateCov,F,Q,lambda,T,stdAngVel)
    % Compute the rotational dynamic model
    angVelEst = prevEstState(10:12);
    [FRot, QRot] = compute_rotational_dynamics(angVelEst, stdAngVel, T);
    
    % Substitute the rotational model
    F(7:12, 7:12) = FRot;
    Q(7:12, 7:12) = QRot;
    
    % Compute predicted state and covariance
    predState = F * prevEstState;
    predQuat = prevEstQuat;
    predStateCov = F * prevEstStateCov * F' + Q;
    
    % Dynamic model for maximum entropy in extent
    predStateCov(13:end, 13:end) = 1/lambda * prevEstStateCov(13:end, 13:end);
end

function [F, Q] = compute_rotational_dynamics(w, std, T)
% Inputs:
%              w:        Angular rate
%              std:     Standard deviation of the angular velocity
%              T:         Sampling time

% Dummy variables
wNorm = norm(w);
if wNorm == 0
    wNorm = 1e-3;
end

S = skew_symmetric_matrix(-w);
c = cos(1/2*T*wNorm);
s = sin(1/2*T*wNorm);
expMat = eye(3,3) + s/wNorm*S + (1-c)/wNorm^2*S^2;

% Construct the state transition matrix
F = [expMat  T*expMat ; zeros(3,3)  eye(3,3)];

% Construct the process noise covariance matrix
G11 = T*eye(3,3) + 2/wNorm^2*(1-c)*S + 1/wNorm^2*(T-2/wNorm*s)*S^2;
G12 = 1/2*T^2*eye(3,3) + 1/wNorm^2*(4/wNorm*s-2*T*c)*S + 1/wNorm^2*(1/2*T^2+2/wNorm*T*s+4/wNorm^2*(c-1))*S^2;
G21 = zeros(3,3);
G22 = T * eye(3,3);
G = [G11 G12; G21 G22];

B = G * [zeros(3,3); eye(3,3)];

% cov = std^2 * diag([0 0 1]);
cov = std^2 * eye(3,3);
Q = B * cov * B';
end

function [out] = skew_symmetric_matrix(x)
out = [0  -x(3)  x(2);...
    x(3)  0  -x(1);...
    -x(2)  x(1)  0];
end
