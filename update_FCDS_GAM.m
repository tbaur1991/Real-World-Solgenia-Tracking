function [X,P,meanInX,meanInY,varInX,varInY] = update_FCDS_GAM(X,P,sig_r,Y,n_upd,nu,nth,artificial_noise,meanInX,meanInY,varInX,varInY,tao)
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
            n_upd = mod(nMeas,n_upd);
        end

        % get parameters of system state
        pos = X(1:3);
        or = X(5);
        R = [cos(or) -sin(or) 0; sin(or) cos(or) 0; 0 0 1];
        h = c1(X(7),0,'lower');
        shape = X(8:end);
    
        % get stacked measurement matrices
        H = zeros(3*n_upd,length(X));
        z_pred = zeros(3*n_upd,1);
        for i = 1:n_upd
            % measurement in local coordinates
            z_loc = R'*(meas(:,i) - pos);
            % get angle parameter of measurement source
            theta = atan2(z_loc(2),z_loc(1));
            % calculate slice of shape
            u = min(max(2*z_loc(3)/h,-1),1);
            % calculate radius
            r = fourier_chebychev_series(shape,theta,u,nu,nth);
            % calculate predicted measurement
            z_pred(3*i-2:3*i) = R*[r*cos(theta);r*sin(theta);u*h/2] + pos - meas(:,i);
            % calculate measurement jacobi matrix
            H(3*i-2:3*i,:) = jacobi_meas_FCDS(theta,u,r,pos,or,h,shape,meas(:,i),nu,nth);
        end
    
        % calculate measurement noise covariance matrix
         if artificial_noise
            % for asymmetric noise
            % first step measurements in local coordinates
            meas_loc = R'*(meas - pos);
            % second step measurements inside or outside
            thetas = atan2(meas_loc(2,:),meas_loc(1,:));
            us = min(max(2*meas_loc(3,:)/h,-1),1);
            rs = zeros(1,n_upd);
            for i = 1:n_upd
                rs(i) = fourier_chebychev_series(shape,thetas(i),us(i),nu,nth);
            end
            inside_outside = vecnorm(meas_loc) - (rs - 3*sig_r);
            % indices of measurements inside
            idx = find(inside_outside < 0);     
            % update estimate of inside measurement mean
            meanInX = 1/(1 + length(idx)/tao)*meanInX + 1/(tao + length(idx))*sum(abs(z_pred(3*idx-2)));
            meanInY = 1/(1 + length(idx)/tao)*meanInY + 1/(tao + length(idx))*sum(abs(z_pred(3*idx-1)));
            z_pred(3*idx-2) = z_pred(3*idx-2) + meanInX; z_pred(3*idx-1) = z_pred(3*idx-1) + meanInY;
            % update estimate of inside measurement variance
            varInX = 1/(1 + length(idx)/tao)*varInX + 1/(tao + length(idx))*sum((z_pred(3*idx-2) - meanInX).^2);
            varInY = 1/(1 + length(idx)/tao)*varInY + 1/(tao + length(idx))*sum((z_pred(3*idx-1) - meanInY).^2);
            % build covariance matrix
            R = sig_r^2*ones(1,3*n_upd);
            R(3*idx-2) = varInX; R(3*idx-1) = varInY; R = diag(R);
        else
            R = sig_r^2*eye(3*n_upd);
        end
    
        % do update
        S = H*P*H' + R;
        K = P*H'/S;
        X = X + K*(zeros(3*n_upd,1) - z_pred);
        P = P - K*H*P;
        P = 0.5*(P+P');
    end
end