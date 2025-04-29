if strcmp(methods,'own')
    % iou elliptic cylinder filter
    iou_cyl(k) = iouEllCyl(X_cyl(:,k),[X_Ref(:,k); h_ref],convRef);
    % iou FCDS filter
    if strcmp(source,'radial')
        if strcmp(filter,'GAM')
            iou_fcds(k) = iouFCDS(X_fcds(:,k),[X_Ref(:,k); h_ref],convRef,nu,nth);
            iou_fcds_full(k) = iouFCDS_full(X_fcds(:,k),[X_Ref(:,k); h_ref],sol_pc,nu,nth,20);
        elseif strcmp(filter,'ERHM')
            X = [X_fcds(1:2,k); X_line(1,k); X_fcds(3:5,k); X_line(2,k); X_fcds(6:end,k)];
            iou_fcds(k) = iouFCDS(X,[X_Ref(:,k); h_ref],convRef,nu,nth);
            iou_fcds_full(k) = iouFCDS_full(X,[X_Ref(:,k); h_ref],sol_pc,nu,nth,20);
        end
    end
    % iou superellipse filter
    iou_super(k) = iouSuper(X_super(:,k),[X_Ref(:,k); h_ref],convRef);
elseif strcmp(methods,'comp')
    % iou SDFS filter
    [iou_sdfs(k),hs_sdfs(k),min_z,max_z] = iouSDFS(X_sdfs(:,k),[X_Ref(:,k); h_ref],convRef,nphi,nth);
    iou_sdfs_full(k) = iouSDFS_full(X_sdfs(:,k),min_z,max_z,[X_Ref(:,k); h_ref],sol_pc,nphi,nth,20);
    % iou SH filter
    [iou_sh(k),hs_sh(k),min_z,max_z] = iouSH(X_sh(:,k),[X_Ref(:,k); h_ref],convRef,nf);
    iou_sh_full(k) = iouSH_full(X_sh(:,k),min_z,max_z,[X_Ref(:,k); h_ref],sol_pc,nf,20);
    % iou 3D GP filter
    [iou_gp(k),hs_gp(k),min_z,max_z] = iouGP(X_gp(:,k),quat_gp(:,k),[X_Ref(:,k); h_ref],convRef,basisAngleArray);
    iou_gp_full(k) = iouGP_full(X_gp(:,k),quat_gp(:,k),min_z,max_z,[X_Ref(:,k); h_ref],sol_pc,basisAngleArray,20);
end