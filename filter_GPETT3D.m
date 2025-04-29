function [ estState, estQuat, estStateCov ] = filter_GPETT3D( predState...
    , predQuat, predStateCov, Y, paramGP, basisAngleArray, eps, n_upd)
% check if number of measurements smaller than n_upd
curNumMeas = size(Y, 1);
numStates = length(predState);
if curNumMeas < n_upd
    n_upd = curNumMeas;
end

% initialize estimated parameters
estState = zeros(size(predState));
estQuat = zeros(size(predQuat));
estStateCov = zeros(size(predStateCov));

% iterate j times over n_upd measurements
for n = 1:ceil(curNumMeas/n_upd)

    % get new measurement set
    if n <= floor(curNumMeas/n_upd)
        measArray = Y(n_upd*n-n_upd+1:n_upd*n,:);
    else
        measArray = Y(end-mod(curNumMeas,n_upd)+1:end,:);
        n_upd = mod(curNumMeas,n_upd);
    end

    % Compute the inverse of P0_extent since it will repeatedly be used in the system dynamic model
    P0_extent = compute_GP_covariance3D(basisAngleArray, basisAngleArray, paramGP);
    P0_extent = P0_extent + eps*eye(size(basisAngleArray,1));       % to prevent numerical errors thrown in matrix inversion
    chol_P0_extent = chol(P0_extent);
    inv_chol_P0_extent = inv(chol_P0_extent);
    inv_P0_extent = inv_chol_P0_extent * inv_chol_P0_extent';
    
    %% Measurement Update
    
    % In the below loop, the following operations are performed for each 3D point measurement:
    % 1. Compute measurement predictions from the predicted state by relying on
    % the nonlinear measurement model.
    % 2. Obtain the corresponding measurement covariance matrix.
    % 3. Linearize the measurement model to employ in EKF equations.
    predMeas = zeros(n_upd * 3, 1);                 	% of the form [x1; y1; z1;...; xn; yn; zn]
    measCov = zeros(n_upd*3);                         	% measurement noise covariance matrix
    linMeasMatrix = zeros(n_upd * 3, numStates);       % linearized measurement model
    for i = 1:n_upd
        iMeas = transpose(measArray(i, :));    % select one measurement
        
        % Obtain predicted measurement and covariance
        [iMeasPredicted,  iMeasCovariance]  = compute_meas_prediction_and_covariance...
            (iMeas, predState, predQuat, paramGP, basisAngleArray, inv_P0_extent);
        
        % Log these variables to later utilize in the update mechanism
        predMeas(1+(i-1)*3 : i*3) = iMeasPredicted;
        measCov(((i-1)*3+1):i*3, ((i-1)*3+1):i*3) = iMeasCovariance;
        
        % Obtain linearized measurement model
        iLinMeasMatrix = compute_linearized_measurement_matrix(iMeas, predState...
            , predQuat, paramGP, basisAngleArray, inv_P0_extent, eps);
        
        % Linearized measurement model for the current measurement
        % linMeasMatrix = [dh/dc, dh/dv = zero, dh/dq, dh/dw = zero, dh/dextent]
        linMeasMatrix(1+(i-1)*3:i*3, :) = iLinMeasMatrix;
    end
    
    % Put the measurements in a column the form: [x1; y1; z1;...; xn; yn; zn]
    tempArray = transpose(measArray);
    curMeasArrayColumn = tempArray(:);
    
    % Realize measurement update
    kalmanGain = predStateCov * linMeasMatrix'...
        / (linMeasMatrix*predStateCov*linMeasMatrix' + measCov);
    estState = predState + kalmanGain * (curMeasArrayColumn - predMeas);
    estStateCov = (eye(numStates) - kalmanGain*linMeasMatrix) * predStateCov;
    estStateCov = (estStateCov + estStateCov')/2;       % to make the covariance matrix symmetric (needed due to numeric errors)
    
    % Quaternion Treatment
    estQuatDev = estState(7:9);
    estQuat = apply_quat_deviation(estQuatDev, predQuat);
    estState(7:9) = zeros(3,1);                         % Reset the orientation deviation

    % rewrite updated state as predicted one
    predState = estState;
    predQuat = estQuat;
    predStateCov = estStateCov;

end

end


function [measPredicted,   measCovariance] = compute_meas_prediction_and_covariance...
    (meas, state, quat, paramGP, basisAngleArray, inv_P0_extent)

meanGP = paramGP{7};    % mean of the GP

% Extract relevant information from the predicted state
center = state(1:3);
extent = state(13:end);

% Compute the following variables to exploit in measurement model
[p, H_f, covMeasBasis, angle_L] = compute_measurement_model_matrices...
    (meas, state, quat, paramGP, basisAngleArray, inv_P0_extent);

% Obtain the measurement prediction by the original nonlinear model
measPredicted = center + p * H_f * extent + p*meanGP*(1-H_f*ones(size(extent,1),1));

% Definition: iCovMeas = k(uk,uk), uk: argument of the current measurement
covGP = compute_GP_covariance3D(angle_L, angle_L, paramGP);
% Definition: iR = k(uk,uk) - k(uk, uf) * inv(K(uf,uf)) * k(uf, uk)
iR_f = covGP - covMeasBasis * inv_P0_extent * covMeasBasis';
% Obtain the covariance of the measurement
stdMeasGP = paramGP{6};
measCovariance = stdMeasGP^2 * eye(3) + p * iR_f * p';
end

function [ linearMeasMat ] = compute_linearized_measurement_matrix( meas, state...
    , quat, paramGP, basisAngleArray, inv_P0_extent, eps)

meanGP = paramGP{7};    % mean of the GP

%% Obtain the measurement prediction by the original nonlinear meas model
[p, H_f, ~, ~] = compute_measurement_model_matrices(meas, state, quat, paramGP...
    , basisAngleArray, inv_P0_extent);
center = state(1:3);
extent = state(13:end);
measPredicted = center + p * H_f * extent + p*meanGP*(1-H_f*ones(size(extent,1),1));

%% Calculate partial derivative wrt target extent
H_extent = p * H_f;     % This partial derivative is taken analytically.

%% Calculate partial derivative wrt linear positions
H_center = zeros(3,3);
for i = 1:3
    stateIncremented = state;   % Initialize the vector
    stateIncremented(i) = stateIncremented(i) + eps;    % increment the state
    centerIncremented = stateIncremented(1:3);          % incremented center
    
    % Compute the output of original measurement model for the incremented state vector
    [p, H_f, ~, ~] = compute_measurement_model_matrices(meas, stateIncremented...
        , quat, paramGP, basisAngleArray, inv_P0_extent);
    measPredictedForIncState = centerIncremented + p * H_f * extent + p*meanGP*(1-H_f*ones(size(extent,1),1));
    
    % Compute numeric derivative
    H_center(:, i) = (measPredictedForIncState - measPredicted) / eps;
end

%% Calculate partial derivative wrt linear velocities
H_linVel = zeros(3,3); % Measurement model does not depend on linear velocities

%% Calculate partial derivative wrt quaternion deviations
H_quatDev = zeros(3,3);
for i = (7:9)
    stateIncremented = state;   % Initialize the vector
    stateIncremented(i) = stateIncremented(i) + eps;
    
    % Compute the output of original measurement model for the incremented state vector
    [p, H_f, ~, ~] = compute_measurement_model_matrices(meas, stateIncremented...
        , quat, paramGP, basisAngleArray, inv_P0_extent);
    measPredictedForIncState = center + p * H_f * extent + p*meanGP*(1-H_f*ones(size(extent,1),1));
    
    % Compute numeric derivative
    H_quatDev(:, i-6) = (measPredictedForIncState - measPredicted) / eps;
end

%% Calculate partial derivative wrt angular velocities
H_angVel = zeros(3,3); % Measurement model does not depend on linear velocities

%% Form the linearized measurement matrix
linearMeasMat = [H_center H_linVel H_quatDev H_angVel H_extent];

end

function [diffUnitVector_G, H_f, covMeasBasis, measAngle_L] = compute_measurement_model_matrices...
    (meas, state, quatPrev, paramGP, basisAngleArray, inv_P0_extent)

% Extract relevant information from the state
center = state(1:3);
quatDev = state(7:9);
quat = apply_quat_deviation(quatDev, quatPrev); % Apply predicted orientation deviation

% Compute the following variables to exploit in measurement model
diffVector_G = meas - center;   % the vector from the center to the measurement
diffVectorMag = norm(diffVector_G);
if diffVectorMag == 0
    diffVectorMag = eps;
end
diffUnitVector_G = diffVector_G / diffVectorMag;

% Express the diffVectorGlobal in Local frame
R_from_G_to_L = rotation_matrix_from_global_to_local(quat);
diffVector_L = R_from_G_to_L * diffUnitVector_G;

% Find the Spherical angle in Local frame corresponding to this vector
[azi_L, elv_L, ~] = cart2sph(diffVector_L(1), diffVector_L(2), diffVector_L(3));
azi_L = mod(azi_L, 2*pi);       % to keep it between 0 and 2*pi
angle_L = [azi_L elv_L];        % angle in the local coordinate frame

% Compute the covariance matrix component from the GP model
covMeasBasis = compute_GP_covariance3D(angle_L, basisAngleArray, paramGP);

H_f = covMeasBasis * inv_P0_extent; % GP model relating the extent and the radial
% function value evaluated at the current measurement angle

measAngle_L = angle_L;  % the angle of the measurement in local coordinate frame

end

function [qOut] = apply_quat_deviation(a, qIn)
qa = 1/ sqrt(4+norm(a)) * [a; 2];   % It depends on the Rodrigues parametrization
qOut = quat_product(qa, qIn);

qOut = qOut/ norm(qOut);    % Normalize due to numeric inprecisions
end

function [out] = quat_product(q1, q2)
q1V = q1(1:3);
q1S = q1(4);
q2V = q2(1:3);
q2S = q2(4);

out = [q1S*q2V + q2S*q1V - cross(q1V,q2V);...
    q1S*q2S - q1V'*q2V];
end
