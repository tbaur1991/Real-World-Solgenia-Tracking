function [X,P] = update_SDFS_UKF(X,P,sig_r,lambda,wm,wc,Y,n_upd,nphi,nth)
    % get dymension of system state
    nx = length(X);

    % check if number of measurements smaller than n_upd
    nMeas = size(Y,2);
    if nMeas < n_upd
        % change number of measurements to be updated and UKF weights
        n_upd = nMeas;
    end
    
    % iterate j times over n_upd measurements
    for n = 1:ceil(nMeas/n_upd)
        % get new measurement set
        if n <= floor(nMeas/n_upd)
            meas = Y(:,n_upd*n-n_upd+1:n_upd*n);
        else
            meas = Y(:,end-mod(nMeas,n_upd)+1:end);
            % change number of measurements to be updated and UKF weights
            n_upd = mod(nMeas,n_upd);
        end

        % calculate sigma points
        A = sqrt(nx + lambda) * chol(P)';
        sig_points = [zeros(size(X)) -A A];
        sig_points = sig_points + repmat(X, 1, size(sig_points, 2));
        numSamples = 2*nx+1;
        
        % predict sigma points for every measurement
        z_predict = zeros(3*n_upd,numSamples);
        for i = 1:numSamples
            % get sample quantities
            pos = sig_points(1:3,i);
            or = sig_points(5,i);
            shape = sig_points(7:end,i);
            % transform measurements to local coordinate system
            R = [cos(or) -sin(or) 0; sin(or) cos(or) 0; 0 0 1];
            z_loc = R'*(meas - pos);
            % calculate orientation vector
            thetas = acos(z_loc(3,:)./vecnorm(z_loc,2)); phis = atan2(z_loc(2,:),z_loc(1,:));
            % calculate radii
            rs = zeros(1,n_upd);
            for j = 1:n_upd
                % calculate radius
                rs(j) = sphDFS_half(thetas(j),phis(j),shape,nphi,nth);
            end
            % calculate predicted measurements
            zp = R*[rs.*sin(thetas).*cos(phis);rs.*sin(thetas).*sin(phis);rs.*cos(thetas)] + pos;
            z_predict(:,i) = zp(:) - meas(:);
        end
        
        % get predicted measurements
        z_pred = sum(wm.*z_predict,2);
        
        % get update matrices
        S = zeros(3*n_upd);
        Psi = zeros(nx,3*n_upd);
        for i = 1:size(sig_points,2)
            S = S + wc(i)*(z_predict(:,i) - z_pred)*(z_predict(:,i) - z_pred)';
            Psi = Psi + wc(i)*(sig_points(1:nx,i) - X)*(z_predict(:,i) - z_pred)';
        end
        S = S + sig_r^2*eye(3*n_upd);
        
        % do UKF measurement update
        K = Psi/S;
        X = X + K*(zeros(3*n_upd,1) - z_pred);
        P = P - K*Psi';
    end
end