function [X,P] = spherical_harmonics_update(X,P,sig_r,Y,wm,wc,lambda,nf,num_coeff,n_upd)
    % get dymension of system state
    nx = length(X);

    % check if number of measurements smaller than n_upd
    nMeas = size(Y,2);
    if nMeas < n_upd
        n_upd = nMeas;
    end
    
    % iterate j times over n_upd measurements
    for n = 1:ceil(nMeas/n_upd)  
        % calculate sigma points
        A = sqrt(nx + lambda) * chol(P)';
        sig_points = [zeros(size(X)) -A A];
        sig_points = sig_points + repmat(X, 1, size(sig_points, 2));
        numSamples = size(sig_points,2);
    
        % get new measurement set
        if n <= floor(nMeas/n_upd)
            meas = Y(:,n_upd*n-n_upd+1:n_upd*n);
        else
            meas = Y(:,end-mod(nMeas,n_upd)+1:end);
            n_upd = mod(nMeas,n_upd);
        end

        % predict measurement for every sample
        z_predict = zeros(3*size(meas,2),numSamples);
        for i = 1:size(meas,2)
            for j = 1:numSamples
                % get parameters of sampled cone
                pos = sig_points(1:3,j);
                or = sig_points(5,j);
                p = sig_points(7:end,j);
                % measurement in local coordinates
                R_rot = [cos(or) -sin(or) 0; sin(or) cos(or) 0; 0 0 1];
                z_loc = R_rot'*(meas(:,i) - pos);
                % calculate orientation vector
                r = norm(z_loc); theta = acos(z_loc(3)/r); phi = atan2(z_loc(2),z_loc(1));
                % calculate unit orientation vector
                pk = (meas(:,i)-pos)/norm(meas(:,i)-pos);
                % calculate radial function value
                radf = sph_harmonics(p,theta,phi,nf,num_coeff);
                % do sample prediction
                z_predict(1+3*(i-1):3*i,j) = pos + pk*radf;
            end
        end
    
        % get predicted measurements
        z_pred = sum(wm.*z_predict,2);
    
        % get update matrices
        R = sig_r^2*eye(3*n_upd);
        S = zeros(3*n_upd);
        Psi = zeros(nx,3*n_upd);
        for i = 1:size(sig_points,2)
            S = S + wc(i)*(z_predict(:,i) - z_pred)*(z_predict(:,i) - z_pred)';
            Psi = Psi + wc(i)*(sig_points(1:nx,i) - X)*(z_predict(:,i) - z_pred)';
        end
        S = S + R;
        
        % do UKF measurement update
        K = Psi/S;
        X = X + K*(meas(:) - z_pred);
        P = P - K*Psi';
    end
end
