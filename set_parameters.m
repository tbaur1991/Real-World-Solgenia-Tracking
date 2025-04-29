% set sensor rate
dt = 0.1;

% set measurement uncertainty
sig_r = 0.025;

% set system uncertainty
sigma_v = 3;
sigma_omega = 5*pi/180;
sigma_ext = 1e-5;

% UKF parameters
alpha_UKF = 0.1; beta_UKF = 2; kappa = 0;

% RAMN parameter
tao = 200;

% covariance transition matrix of process noise
G = zeros(6,2);
G(1,1) = dt^2/2; G(2,1) = dt^2/2; G(3,1) = dt^2/2; 
G(4,1) = dt; G(5,2) = dt^2/2; G(6,2) = dt; 

% get height of reference boat
h_ref = max(sol_mesh.Vertices(:,3))*2;

% get 2D convex hull of boat reference
pc = pcdownsample(pointCloud(sol_pc'),'gridAverage',0.1);
dsPoints = pc.Location;
c = convhull(dsPoints(:,1),dsPoints(:,2));
convRef = dsPoints(c,1:2)';

% filter parameters
if strcmp(methods,'own')
    %% Elliptic cylinder parameters
    % memory allocation
    nx = 9; % [x y z v or omega a b h]
    X_cyl = zeros(nx,nSamples);
    P_cyl = zeros(nx,nx,nSamples);
    iou_cyl = zeros(1,nSamples);
    
    % sigma point parameters
    if strcmp(filter,'ERHM')
        nxu = nx + n_upd;
    elseif strcmp(filter,'GAM')
        nxu = nx;
    else
        error('Wrong filter setting!')
    end
    lambda_cyl = alpha_UKF^2*(nxu + kappa) - nxu;
    
    % calculate weights for sigma points
    wm_cyl(1) = lambda_cyl/(nxu + lambda_cyl);
    wc_cyl(1) = lambda_cyl/(nxu + lambda_cyl) + (1 - alpha_UKF^2 + beta_UKF);
    wm_cyl(2:2*nxu + 1) = 1/(2*(nxu + lambda_cyl));
    wc_cyl(2:2*nxu + 1) = 1/(2*(nxu + lambda_cyl));

    % elliptic cylinder RAMN parameter
    meanX_cyl = 0; meanY_cyl = 0;
    varX_cyl = 0; varY_cyl = 0;

    % elliptic cylinder covariance transition matrix of process noise
    G_cyl = blkdiag(G,dt*eye(3));

    %% FCDS parameters
    % set number of coefficients for FCDS shape expansion
    nu = 5;
    nth = 7;
    nCoeff = 1 + nu + nth + nu*nth;
    
    % System state and evaluation allocation for cylindric shape estimation
    if strcmp(filter,'GAM')
        nx = 7 + nCoeff; % [x y z v or omega h shape]
    elseif strcmp(filter,'ERHM')
        nx = 5 + nCoeff; % [x y v or omega shape]
        % line system state
        X_line = zeros(2,nSamples);
        P_line = zeros(2,2,nSamples);
        % sigma point parameters line
        nx_line = 2 + n_upd;
        lambda_line = alpha_UKF^2*(nx_line + kappa) - nx_line;
        % calculate weights for sigma points
        wm_line(1) = lambda_line/(nx_line + lambda_line);
        wc_line(1) = lambda_line/(nx_line + lambda_line) + (1 - alpha_UKF^2 + beta_UKF);
        wm_line(2:2*nx_line + 1) = 1/(2*(nx_line + lambda_line));
        wc_line(2:2*nx_line + 1) = 1/(2*(nx_line + lambda_line));
    end
    X_fcds = zeros(nx,nSamples);
    P_fcds = zeros(nx,nx,nSamples);
    iou_fcds = zeros(1,nSamples);
    iou_fcds_full = zeros(1,nSamples);

    % FCDS RAMN parameter
    meanX_fcds = 0; meanY_fcds = 0;
    varX_fcds = 0; varY_fcds = 0;

    % FCDS covariance transition matrix of process noise
    G_fcds = blkdiag(G,dt*eye(nCoeff+1));

    %% Superellipse parameters
    % System state and measurement memory allocation cylindric shape estimation
    nx = 11; % [x y z v or omega a b h epsilon ty]
    X_super = zeros(nx,nSamples);
    P_super = zeros(nx,nx,nSamples);
    iou_super = zeros(1,nSamples);
    
    % sigma point parameters
    if strcmp(filter,'ERHM')
        nxu = nx + n_upd;
    elseif strcmp(filter,'GAM')
        nxu = nx;
    end
    lambda_super = alpha_UKF^2*(nxu + kappa) - nxu;
    
    % calculate weights for sigma points
    wm_super(1) = lambda_super/(nxu + lambda_super);
    wc_super(1) = lambda_super/(nxu + lambda_super) + (1 - alpha_UKF^2 + beta_UKF);
    wm_super(2:2*nxu + 1) = 1/(2*(nxu + lambda_super));
    wc_super(2:2*nxu + 1) = 1/(2*(nxu + lambda_super));

    % superellipse RAMN parameter
    meanX_super = 0; meanY_super = 0;
    varX_super = 0; varY_super = 0;

    % superellipse covariance transition matrix of process noise
    G_super = blkdiag(G,dt*eye(5));
elseif strcmp(methods,'comp')
    %% SDFS parameters
    % set number of coefficients for SDFS shape expansion
    nphi = 5;
    nth = 5;
    nCoeff_sdfs = 1 + nth + nth*nphi;
    
    % System state and evaluation allocation for cylindric shape estimation
    nx = 6 + nCoeff_sdfs;
    X_sdfs = zeros(nx,nSamples);
    P_sdfs = zeros(nx,nx,nSamples);

    % evaluation memory allocation
    iou_sdfs = zeros(1,nSamples);
    iou_sdfs_full = zeros(1,nSamples);
    hs_sdfs = zeros(1,nSamples);
    
    % sigma point parameters FCDS
    lambda_sdfs = alpha_UKF^2*(nx + kappa) - nx;
    
    % calculate weights for sigma points
    wm_sdfs(1) = lambda_sdfs/(nx + lambda_sdfs);
    wc_sdfs(1) = lambda_sdfs/(nx + lambda_sdfs) + (1 - alpha_UKF^2 + beta_UKF);
    wm_sdfs(2:2*nx + 1) = 1/(2*(nx + lambda_sdfs));
    wc_sdfs(2:2*nx + 1) = 1/(2*(nx + lambda_sdfs));

    % SDFS covariance transition matrix of process noise
    G_sdfs = blkdiag(G,dt*eye(nCoeff_sdfs));

    %% SH parameters
    % set number of coefficients for shape expansion
    nf = 5;
    nCoeff_sh = (2*nf + 1)^2;
    
    % System state and measurement memory allocation cylindric shape tracking
    nx_sh = 6 + nCoeff_sh;
    X_sh = zeros(nx_sh,nSamples);
    P_sh = zeros(nx_sh,nx_sh,nSamples);
    
    % evaluation memory allocation
    iou_sh = zeros(1,nSamples);
    iou_sh_full = zeros(1,nSamples);
    hs_sh = zeros(1,nSamples);
    
    % calculate weights for sigma points of spherical harmonics
    lambda_sh = alpha_UKF^2*(nx_sh + kappa) - nx_sh;
    wm_sh(1) = lambda_sh/(nx_sh + lambda_sh);
    wc_sh(1) = lambda_sh/(nx_sh + lambda_sh) + (1 - alpha_UKF^2 + beta_UKF);
    wm_sh(2:2*nx_sh + 1) = 1/(2*(nx_sh + lambda_sh));
    wc_sh(2:2*nx_sh + 1) = 1/(2*(nx_sh + lambda_sh));

    % SH covariance transition matrix of process noise
    G_sh = blkdiag(G,dt*eye(nCoeff_sh));

    %% GP parameters
    % Gaussian Process (GP) Parameters
    meanGP = 0;                     % mean of the GP
    stdPriorGP = 1;                 % prior standard deviation of the GP
    stdRadiusGP = 0.2;              % standard deviation of the radius
    scaleLengthGP = pi/14;           % lengthscale
    kernelTypeGP = 1;               % 1:Squared exponential , 2: Matern3_2, 3: Matern5_2
    isObjSymmetric = 0;             % it is used to adjust the kernel function accordingly
    stdMeasGP = sig_r;                % standard deviation of the measurements (used in the GP model)
    paramGP = {kernelTypeGP, stdPriorGP, stdRadiusGP, scaleLengthGP, isObjSymmetric, stdMeasGP, meanGP};
    
    % Determine the Basis Angles using Fibonacci lattice
    n = 200;
    lats = zeros(2*n+1,1);
    lons = zeros(2*n+1,1);
    gol_rat = (1+sqrt(5))/2;
    for i = -n:n
        lats(i+n+1) = asin((2*i)/(2*n+1))*180/pi;
        lons(i+n+1) = mod(i,gol_rat)*360/gol_rat;
        if lons(i+n+1) < -180
            lons(i+n+1) = 360 + lons(i+n+1);
        elseif lons(i+n+1) > 180
            lons(i+n+1) = lons(i+n+1) - 360;
        end
    end
    lats = lats*pi/180; 
    lons = lons*pi/180; 
    basisAngleArray = [lons lats];
    basisAngleArray = sortrows(basisAngleArray, 2);
    basisAngleArray = sortrows(basisAngleArray, 1);
    numBasisAngles = size(basisAngleArray,1);
    
    % System state and measurement memory allocation cylindric shape tracking
    nx = 12 + numBasisAngles;
    X_gp = zeros(nx,nSamples);
    quat_gp = zeros(4,nSamples);
    P_gp = zeros(nx,nx,nSamples);
    
    % EKF Paramaters
    stdCenter = 1e-1;               % std dev of the process noise for object center
    stdAngVel = 1e-1;               % std dev of the process noise for angular velocity
    lambda_gp = 0.99;  

    % define the Process Model
    F_lin = kron([1 dt; 0 1], eye(3));       % constant velocity model for linear motion
    F_rot = eye(6,6);                       % a substitute matrix for the rotational dynamics (it will be computed in the filter at each recursion)
    F_f = eye(numBasisAngles);      
    F = blkdiag(F_lin, F_rot, F_f);
    
    Q_lin = kron([dt^3/3 dt^2/2; dt^2/2 dt], stdCenter^2*eye(3));   % process covariance of linear motion
    Q_rot = eye(6,6);                                           % a substitute matrix for the covariance of the rotational motion (it will be computed in the filter at each recursion)
    Q_extent = zeros(numBasisAngles);                           % predicted covariance of the extent will be computed within the filter according to Eqn. 20 
    Q = blkdiag(Q_lin, Q_rot, Q_extent);
    
    % evaluation memory allocation
    iou_gp = zeros(1,nSamples);
    iou_gp_full = zeros(1,nSamples);
    hs_gp = zeros(1,nSamples);
else
    error('Wrong methods setting!')
end