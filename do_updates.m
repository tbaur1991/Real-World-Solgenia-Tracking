if strcmp(methods,'own')
    % update elliptic cylinder filter
    [X_cyl(:,k),P_cyl(:,:,k),meanX_cyl,meanY_cyl,varX_cyl,varY_cyl] = elliptic_cylinder_update(X_cyl(:,k),P_cyl(:,:,k),meas{k,1}(:,idxp),...
     sig_r,lambda_cyl,wm_cyl,wc_cyl,n_upd,meanX_cyl,meanY_cyl,varX_cyl,varY_cyl,tao,artificial_noise,filter,source,alpha_UKF,beta_UKF,kappa);
    % update FCDS filter
    if strcmp(source,'radial')
        if strcmp(filter,'GAM')
            [X_fcds(:,k),P_fcds(:,:,k),meanX_fcds,meanY_fcds,varX_fcds,varY_fcds] = update_FCDS_GAM(X_fcds(:,k),P_fcds(:,:,k),sig_r,...
                                         meas{k,1}(:,idxp),n_upd,nu,nth,artificial_noise,meanX_fcds,meanY_fcds,varX_fcds,varY_fcds,tao);
        elseif strcmp(filter,'ERHM')
            [X_line(:,k),P_line(:,:,k)] = update_line_ERHM(X_line(:,k),P_line(:,:,k),sig_r,double(meas{k,1}(:,idxp)),n_upd,lambda_line,wm_line,wc_line);
            [X_fcds(:,k),P_fcds(:,:,k),meanX_fcds,meanY_fcds,varX_fcds,varY_fcds] = update_FCDS_ERHM(X_fcds(:,k),P_fcds(:,:,k),X_line(:,k),...
                                          sig_r,meas{k,1}(:,idxp),n_upd,nu,nth,artificial_noise,meanX_fcds,meanY_fcds,varX_fcds,varY_fcds,tao);
        end
    end
    % update superellipse filter
    [X_super(:,k),P_super(:,:,k),meanX_super,meanY_super,varX_super,varY_super] = superellipse3D_UKF_update(X_super(:,k),P_super(:,:,k),meas{k,1}(:,idxp),...
        sig_r,lambda_super,wm_super,wc_super,n_upd,meanX_super,meanY_super,varX_super,varY_super,tao,artificial_noise,filter,source,alpha_UKF,beta_UKF,kappa);
elseif strcmp(methods,'comp')
    % update SDFS filter
    [X_sdfs(:,k),P_sdfs(:,:,k)] = update_SDFS_UKF(X_sdfs(:,k),P_sdfs(:,:,k),sig_r,lambda_sdfs,wm_sdfs,wc_sdfs,meas{k,1}(:,idxp),n_upd,nphi,nth);
    % update SH filter
    [X_sh(:,k),P_sh(:,:,k)] = spherical_harmonics_update(X_sh(:,k),P_sh(:,:,k),sig_r,meas{k,1}(:,idxp),wm_sh,wc_sh,lambda_sh,nf,nCoeff_sh,n_upd);
    % update 3D GP filter
    [X_gp(:,k),quat_gp(:,k),P_gp(:,:,k)] = filter_GPETT3D(X_gp(:,k),quat_gp(:,k),P_gp(:,:,k),meas{k,1}(:,idxp)',paramGP,basisAngleArray,1e-6,n_upd);
end