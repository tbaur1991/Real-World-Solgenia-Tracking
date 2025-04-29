if strcmp(methods,'own')
    %% elliptic cylinder filter plot
    % prepare 3D plot
    figure(1)
    tiledlayout(3,1)
    nexttile
    ax = gca;
    set(gcf,'Color','w')
    % plot measurements
    plot3(meas{k,1}(1,:),meas{k,1}(2,:),meas{k,1}(3,:),'.r','MarkerSize',8);
    axis equal; hold on
    % plot solgenia reference
    pose_ref = X_Ref(:,k);
    ref = rotate(sol_mesh,[pose_ref(5)*180/pi 0 0]);
    ref = translate(ref,(pose_ref(1:3))');
    show(ref)
    % plot estimate
    syms theta u
    pos = X_cyl(1:3,k); or = X_cyl(5,k); 
    a = c1(X_cyl(7,k),0,'lower'); b = c1(X_cyl(8,k),0,'lower'); h = c1(X_cyl(9,k),0,'lower');
    x = a*cos(theta);
    y = b*sin(theta);
    z = u*h/2;
    p2 = fsurf(x,y,z,[0 2*pi -1 1],'FaceAlpha',0.6,'FaceColor','#4DBEEE','MeshDensity',15);
    % translate surf plot
    t = hgtransform('Parent',ax); set(p2,'Parent',t);
    m = makehgtform('translate',pos); set(t,'Matrix',m)
    % rotate surf plot
    t = hgtransform('Parent',t); set(p2,'Parent',t);
    Rz = makehgtform('zrotate',or); set(t,'Matrix',Rz);
    % axes
    set(gca,'TickLength',[0 0])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    xlabel('$x/m$','Interpreter','latex')
    ylabel('$y/m$','Interpreter','latex')
    zlabel('$z/m$','Interpreter','latex')
    % view of figure
    view(mod(pose_ref(5)*180/pi,360)-25,30)
    % draw
    drawnow
    hold off

    %% FCDS filter plot
    % prepare 3D plot
    nexttile
    ax = gca;
    % plot measurements
    plot3(meas{k,1}(1,:),meas{k,1}(2,:),meas{k,1}(3,:),'.r','MarkerSize',8);
    axis equal; hold on
    % plot solgenia reference
    pose_ref = X_Ref(:,k);
    ref = rotate(sol_mesh,[pose_ref(5)*180/pi 0 0]);
    ref = translate(ref,(pose_ref(1:3))');
    show(ref)
    % plot estimate 
    syms theta u
    if strcmp(filter,'GAM')
        pos = X_fcds(1:3,k); or = X_fcds(5,k); 
        h = c1(X_fcds(7,k),0,'lower');
        shape = X_fcds(8:end,k);
    elseif strcmp(filter,'ERHM')
        pos = [X_fcds(1:2,k); X_line(1,k)]; or = X_fcds(4,k); 
        h = c1(X_line(2,k),0,'lower');
        shape = X_fcds(6:end,k);
    end
    r = fourier_chebychev_series(shape,theta,u,nu,nth);
    x = r*cos(theta);
    y = r*sin(theta);
    z = u*h/2;
    p2 = fsurf(x,y,z,[0 2*pi -1 1],'FaceAlpha',0.6,'FaceColor','#4DBEEE','MeshDensity',15);
    % translate surf plot
    t = hgtransform('Parent',ax); set(p2,'Parent',t);
    m = makehgtform('translate',pos); set(t,'Matrix',m)
    % rotate surf plot
    t = hgtransform('Parent',t); set(p2,'Parent',t);
    Rz = makehgtform('zrotate',or); set(t,'Matrix',Rz);
    % axes
    set(gca,'TickLength',[0 0])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    xlabel('$x/m$','Interpreter','latex')
    ylabel('$y/m$','Interpreter','latex')
    zlabel('$z/m$','Interpreter','latex')
    % view of figure
    view(mod(pose_ref(5)*180/pi,360)-25,30)
    % draw
    drawnow
    hold off

    %% superellipse filter plot
    % prepare 3D plot
    nexttile
    ax = gca;
    set(gcf,'Color','w')
    % plot measurements
    p1 = plot3(meas{k,1}(1,:),meas{k,1}(2,:),meas{k,1}(3,:),'.r','MarkerSize',8);
    axis equal; hold on
    % plot solgenia reference
    pose_ref = X_Ref(:,k);
    ref = rotate(sol_mesh,[pose_ref(5)*180/pi 0 0]);
    ref = translate(ref,(pose_ref(1:3))');
    show(ref)
    % plot estimate
    syms theta u
    pos = X_super(1:3,k); or = X_super(5,k); 
    a = c1(X_super(7,k),0,'lower'); b = c1(X_super(8,k),0,'lower'); h = c1(X_super(9,k),0,'lower');
    e = c1(X_super(10,k),1,'lower'); ty = c2(X_super(11,k),-1,1);
    x = a*sign(cos(theta)).*abs(cos(theta)).^(2/e);
    y = b*sign(sin(theta)).*abs(sin(theta)).^(2/e);
    z = u*h/2;
    % tapering y
    y = (ty*x/a + 1).*y;
    p2 = fsurf(x,y,z,[0 2*pi -1 1],'FaceAlpha',0.6,'FaceColor','#4DBEEE','MeshDensity',15);
    % translate surf plot
    t = hgtransform('Parent',ax); set(p2,'Parent',t);
    m = makehgtform('translate',pos); set(t,'Matrix',m)
    % rotate surf plot
    t = hgtransform('Parent',t); set(p2,'Parent',t);
    Rz = makehgtform('zrotate',or); set(t,'Matrix',Rz);
    % axes
    set(gca,'TickLength',[0 0])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    xlabel('$x/m$','Interpreter','latex')
    ylabel('$y/m$','Interpreter','latex')
    zlabel('$z/m$','Interpreter','latex')
    % legend
    leg = legend([p1,p2],{'Measurement','Estimate'});
    leg.Interpreter = "latex";
    leg.Orientation = "horizontal";
    leg.NumColumns = 3;
    leg.Box = "off";
    leg.Layout.Tile = "south";
    leg.FontSize = 12;
    % view of figure
    view(mod(pose_ref(5)*180/pi,360)-25,30)
    % draw
    drawnow
    hold off
    % set position
    set(gcf,'Position',[2520,67,631,911])
elseif strcmp(methods,'comp')
    %% SDFS filter plot
    % prepare 3D plot
    figure(1)
    tiledlayout(1,1)
    nexttile
    ax = gca;
    set(gcf,'Color','w')
    % plot measurements
    p1 = plot3(meas{k,1}(1,:),meas{k,1}(2,:),meas{k,1}(3,:),'.r','MarkerSize',8);
    axis equal; hold on
    % plot solgenia reference
    pose_ref = X_Ref(:,k);
    ref = rotate(sol_mesh,[pose_ref(5)*180/pi 0 0]);
    ref = translate(ref,(pose_ref(1:3))');
    show(ref)
    % plot estimate
    syms theta phi 
    pos = X_sdfs(1:3,k);
    or = X_sdfs(5,k);
    p = X_sdfs(7:end,k);
    r = sphDFS_half(theta,phi,p,nphi,nth);
    x = r*sin(theta)*cos(phi);
    y = r*sin(theta)*sin(phi);
    z = r*cos(theta);
    p2 = fsurf(x,y,z,[0 2*pi 0 pi],'FaceAlpha',0.6,'FaceColor','#4DBEEE');
    % translate surf plot
    t = hgtransform('Parent',ax);
    set(p2,'Parent',t);
    m = makehgtform('translate',pos);
    set(t,'Matrix',m)
    % rotate surf plot
    t = hgtransform('Parent',t);
    set(p2,'Parent',t);
    Rz = makehgtform('zrotate',or);
    set(t,'Matrix',Rz);
    % axes
    set(gca,'TickLength',[0 0])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('$z$','Interpreter','latex')
    % legend
    leg = legend([p1,p2],{'Measurement','Estimate'});
    leg.Interpreter = "latex";
    leg.Orientation = "horizontal";
    leg.NumColumns = 3;
    leg.Box = "off";
    leg.Layout.Tile = "south";
    leg.FontSize = 12;
    % view of figure
    view(mod(pose_ref(5)*180/pi,360)-25,30)
    % draw
    drawnow
    hold off

    %% SH filter plot
    % prepare 3D plot
    figure(2)
    tiledlayout(1,1)
    nexttile
    ax = gca;
    set(gcf,'Color','w')
    % plot measurements
    p1 = plot3(meas{k,1}(1,:),meas{k,1}(2,:),meas{k,1}(3,:),'.r','MarkerSize',8);
    axis equal; hold on
    % plot solgenia reference
    pose_ref = X_Ref(:,k);
    ref = rotate(sol_mesh,[pose_ref(5)*180/pi 0 0]);
    ref = translate(ref,(pose_ref(1:3))');
    show(ref)
    % plot estimate
    syms theta phi 
    pos = X_sh(1:3,k);
    or = X_sh(5,k);
    p = X_sh(7:end,k);
    r = sph_harmonics(p,theta,phi,nf,nCoeff_sh);
    x = r*sin(theta)*cos(phi);
    y = r*sin(theta)*sin(phi);
    z = r*cos(theta);
    p2 = fsurf(x,y,z,[0 2*pi 0 pi],'FaceAlpha',0.6,'FaceColor','#4DBEEE');
    % translate surf plot
    t = hgtransform('Parent',ax);
    set(p2,'Parent',t);
    m = makehgtform('translate',pos);
    set(t,'Matrix',m)
    % rotate surf plot
    t = hgtransform('Parent',t);
    set(p2,'Parent',t);
    Rz = makehgtform('zrotate',or);
    set(t,'Matrix',Rz);
    % axes
    set(gca,'TickLength',[0 0])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('$z$','Interpreter','latex')
    % legend
    leg = legend([p1,p2],{'Measurement','Estimate'});
    leg.Interpreter = "latex";
    leg.Orientation = "horizontal";
    leg.NumColumns = 3;
    leg.Box = "off";
    leg.Layout.Tile = "south";
    leg.FontSize = 12;
    % view of figure
    view(mod(pose_ref(5)*180/pi,360)-25,30)
    % draw
    drawnow
    hold off

    %% 3D GP filter plot
    % prepare 3D plot
    figure(3)
    tiledlayout(1,1)
    nexttile
    ax = gca;
    set(gcf,'Color','w')
    % plot measurements
    p1 = plot3(meas{k,1}(1,:),meas{k,1}(2,:),meas{k,1}(3,:),'.r','MarkerSize',8);
    axis equal; hold on
    % plot solgenia reference
    pose_ref = X_Ref(:,k);
    ref = rotate(sol_mesh,[pose_ref(5)*180/pi 0 0]);
    ref = translate(ref,(pose_ref(1:3))');
    show(ref)
    % plot estimate
    pos = X_gp(1:3,k);
    quat = quat_gp(:,k);
    shape = X_gp(13:end,k);
    % Calculate the rotation matrix from Local to Global by the quternions
    R = rotation_matrix_from_global_to_local(quat); R = R';
    % calculate points in cartesian coordinates
    [x,y,z] = sph2cart(basisAngleArray(:,1),basisAngleArray(:,2),shape);
    pts = [x';y';z'];
    % Transform the extent to the Global frame
    pts = pos + R*pts;
    % plot surface
    T = boundary(pts', 0.5);
    p2 = trisurf(T,pts(1,:)',pts(2,:)',pts(3,:)','FaceAlpha',0.6,'FaceColor','#4DBEEE');
    % axes
    set(gca,'TickLength',[0 0])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    zlabel('$z$','Interpreter','latex')
    % legend
    leg = legend([p1,p2],{'Measurement','Estimate'});
    leg.Interpreter = "latex";
    leg.Orientation = "horizontal";
    leg.NumColumns = 3;
    leg.Box = "off";
    leg.Layout.Tile = "south";
    leg.FontSize = 12;
    % view of figure
    view(mod(pose_ref(5)*180/pi,360)-25,30)
    % draw
    drawnow
    hold off
end