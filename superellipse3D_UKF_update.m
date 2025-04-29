function [X,P,meanInX,meanInY,varInX,varInY] = superellipse3D_UKF_update(X,P,Y,sig_r,lambda,wm,wc,n_upd,meanInX,meanInY,varInX,varInY,tao,artificial_noise,filter,source,alpha,beta,kappa)
    % get dymension of system state
    nx = length(X); 

    % check if number of measurements smaller than n_upd
    nMeas = size(Y,2);
    if nMeas < n_upd
        % change number of measurements to be updated and UKF weights
        n_upd = nMeas;
        if strcmp(filter,'ERHM')
            nxu = nx + n_upd;
        elseif strcmp(filter,'GAM')
            nxu = nx;
        else
            nxu = 0;
        end
        lambda = alpha^2*(nxu + kappa) - nxu;
        wm = zeros(1,2*nxu+1); wc = wm;
        wm(1) = lambda/(nxu + lambda);
        wc(1) = lambda/(nxu + lambda) + (1 - alpha^2 + beta);
        wm(2:2*nxu + 1) = 1/(2*(nxu + lambda));
        wc(2:2*nxu + 1) = 1/(2*(nxu + lambda));
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
            if strcmp(filter,'ERHM')
                nxu = nx + n_upd;
            elseif strcmp(filter,'GAM')
                nxu = nx;
            else
                nxu = 0;
            end
            lambda = alpha^2*(nxu + kappa) - nxu;
            wm = zeros(1,2*nxu+1); wc = wm;
            wm(1) = lambda/(nxu + lambda);
            wc(1) = lambda/(nxu + lambda) + (1 - alpha^2 + beta);
            wm(2:2*nxu + 1) = 1/(2*(nxu + lambda));
            wc(2:2*nxu + 1) = 1/(2*(nxu + lambda));
        end

        % calculate sigma points
        if strcmp(filter,'ERHM')
            Xu = [X; zeros(n_upd,1)];
            L = blkdiag(chol(P)',eye(n_upd));
            nxu = nx + n_upd;
        elseif strcmp(filter,'GAM')
            Xu = X;
            L = chol(P)';
            nxu = nx;
        else
            error('Wrong filter setting!')
        end
        A = sqrt(nxu + lambda) * L;
        sig_points = [zeros(size(Xu)) -A A];
        sig_points = sig_points + repmat(Xu, 1, size(sig_points, 2));
        sig_points(nx+1:end,:) = 0.5*(1 + erf(sig_points(nx+1:end,:)/sqrt(2)));
        numSamples = size(sig_points,2);
        
        % predict sigma points for every measurement
        z_predict = zeros(3*n_upd,numSamples);
        for i = 1:numSamples
            % get sample quantities
            pos = sig_points(1:3,i);
            or = sig_points(5,i);
            a = c1(sig_points(7,i),0,'lower'); 
            b = c1(sig_points(8,i),0,'lower'); 
            h = c1(sig_points(9,i),0,'lower');
            e = c1(sig_points(10,i),1,'lower'); 
            ty = c2(sig_points(11,i),-1,1);
            % transform measurements to local coordinate system
            R = [cos(or) -sin(or) 0; sin(or) cos(or) 0; 0 0 1];
            meas_loc = R'*(meas - pos);
            meas_loc(3,:) = abs(meas_loc(3,:));
            % do inverse tapering transformation
            meas_loc(2,:) = meas_loc(2,:)./(ty*meas_loc(1,:)./a + 1);
            % extrusion factors
            if strcmp(filter,'ERHM')
                us = sig_points(12:end,i);
            elseif strcmp(filter,'GAM')
                us = min(max(2*meas_loc(3,:)/h,0),1)';
            else
                us = 0;
            end
            % do sample prediction with ray function or polygonal chain
            % approximation
            if strcmp(source,'radial')
                % calculate intersection point
                xi = sign(meas_loc(1,:)).*power(complex(1./(1/(abs(a)^e) + abs(power(complex((meas_loc(2,:)./(meas_loc(1,:)*b))),e)))),1/e);
                yi = xi.*(meas_loc(2,:)./meas_loc(1,:));
                zs = [xi; yi];
            elseif strcmp(source,'projected')
                % calculate polygonal chain and projection point
                ps = as_polygon_super(sig_points(:,i),25);
                zs = project(ps, meas_loc(1:2,:)); 
            else
                error('Wrong source association method chosen.')
            end   
            zp = [zs; us'*h/2];
            % measurement prediction
            z_predict(:,i) = real(zp(:) - meas_loc(:));
        end
        
        % get predicted measurements
        z_pred = sum(wm.*z_predict,2);
    
        % calculate measurement noise covariance matrix
        if artificial_noise
            % for asymmetric noise
            % first step measurements in local coordinates
            pos = X(1:3); or = X(5); 
            R = [cos(or) -sin(or) 0; sin(or) cos(or) 0; 0 0 1]; 
            a = c1(X(7),0,'lower'); b = c1(X(8),0,'lower'); 
            e = c1(X(10),1,'lower'); ty = c2(X(11),-1,1);
            meas_loc = R'*(meas - pos);
            meas_loc(2,:) = meas_loc(2,:)./(ty*meas_loc(1,:)./a + 1);
            % second step measurements inside or outside
            inside_outside = abs(meas_loc(1,:)/(a - 3*sig_r)).^e + abs(meas_loc(2,:)/(b - 3*sig_r)).^e - 1;
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
        
        % get update matrices
        S = zeros(3*n_upd);
        Psi = zeros(nx,3*n_upd);
        for i = 1:size(sig_points,2)
            S = S + wc(i)*(z_predict(:,i) - z_pred)*(z_predict(:,i) - z_pred)';
            Psi = Psi + wc(i)*(sig_points(1:nx,i) - X)*(z_predict(:,i) - z_pred)';
        end
        S = S + R;
        
        % do UKF measurement update
        K = Psi/S;
        X = X + K*(zeros(3*n_upd,1) - z_pred);
        P = P - K*Psi';
    end
end