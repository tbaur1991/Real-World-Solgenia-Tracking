%% get parameters for state initialization
% init position
pos = zeros(3,1);
pos(1) = min(meas{k,1}(1,:)) + (max(meas{k,1}(1,:)) - min(meas{k,1}(1,:)))/2;
pos(2) = min(meas{k,1}(2,:)) + (max(meas{k,1}(2,:)) - min(meas{k,1}(2,:)))/2;
pos(3) = min(meas{k,1}(3,:));

% init orientation
c = cov(meas{k,1}(1,:),meas{k,1}(2,:));
[v,d] = eig(c); [~,i] = max(diag(d));
or = atan2(v(2,i),v(1,i));

% init shape
R_rot = [cos(or) -sin(or) 0; sin(or) cos(or) 0; 0 0 1];
mm = R_rot'*(meas{k,1} - pos);
a = max(abs(mm(1,:)));
b = abs(max(mm(2,:)));
h = max(meas{k,1}(3,:)) - min(meas{k,1}(3,:));

%% define initial states
if strcmp(methods,'own')
    % init elliptic cylinder filter
    X_cyl(:,k) = [pos(1:2); pos(3)+h/2; 0; or; 0; a; b; h];
    P_cyl(:,:,k) = blkdiag(c, 2, 2, 10*pi/180*eye(2), 2*eye(3)); 
    % init FCDS filter
    if strcmp(filter,'GAM')
        X_fcds(:,k) = [pos(1:2); pos(3)+h/2; 0; or; 0; h; 4*b; zeros(nCoeff-1,1)];
        P_fcds(:,:,k) = blkdiag(c, 2, 2, 10*pi/180*eye(2), 2*eye(nCoeff+1));
    elseif strcmp(filter,'ERHM')
        X_fcds(:,k) = [pos(1:2); 0; or; 0; 4*b; zeros(nCoeff-1,1)];
        P_fcds(:,:,k) = blkdiag(c, 2, 10*pi/180*eye(2), 2*eye(nCoeff));
        X_line(:,k) = [pos(3)+h/2; h];
        P_line(:,:,k) = blkdiag(2, 0.5);
    end
    % init superellipse filter
    X_super(:,k) = [pos(1:2); pos(3)+h/2; 0; or; 0; a; b; h; 3; 0];
    P_super(:,:,k) = blkdiag(c, 2, 2, 10*pi/180*eye(2), 2*eye(3), 1, 1/3); 
elseif strcmp(methods,'comp')
    % init SDFS filter
    X_sdfs(:,k) = [pos(1:2);pos(3)+h/2;0;or;0;b;zeros(nCoeff_sdfs-1,1)];
    P_sdfs(:,:,k) = blkdiag(c,.2*eye(2),5*pi/180*eye(2),0.3*eye(nCoeff_sdfs));
    % init SH filter
    X_sh(:,k) = [pos(1:2);pos(3)+h/2;0;or;0;b;zeros(nCoeff_sh-1,1)];
    P_sh(:,:,k) = blkdiag(c,.2*eye(2),5*pi/180*eye(2),0.3*eye(nCoeff_sh));
    % init 3D GP filter
    % initial mean
    c0 = pos;        % initial linear velocity
    v0 = eps*ones(3,1);                          % initial linear velocity
    a0 = [0 0 0]';                                                  % initial orientation deviation
    w0 = [0 0 0]';                                                  % initial angular velocity
    f0 = meanGP * ones(numBasisAngles, 1);                          % initial extent (is initialized regarding the GP model)
    X_gp(:,k) = [c0; v0; a0; w0; f0];
    quat_gp(:,k) = [0 0 0 1]';                                      % initial quaternions
    % Initial covariance
    P0_c = blkdiag(c,0.2);
    P0_v = 0.2 * eye(3);
    P0_a = 1e-5 * eye(3,3);
    P0_w = 1e-0 * eye(3,3);
    P0_extent = compute_GP_covariance3D(basisAngleArray, basisAngleArray, paramGP);  
    P_gp(:,:,k) = blkdiag(P0_c, P0_v, P0_a, P0_w, P0_extent);
end
