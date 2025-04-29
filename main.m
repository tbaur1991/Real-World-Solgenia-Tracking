% *************************************************************************
%    I   SSSS   DDDD
%    I   S      D   D 
%    I   SSSS   D   D    Institut fr Systemdynamik
%    I      S   D   D
%    I   SSSS   DDDD   
% *************************************************************************
% Author: Tim Baur
% University: University of Applied Sciences "HTWG Konstanz"
% E-Mail: tbaur@htwg-konstanz.de    
% Creation date: 29.04.2025
% File name: main.m
% *************************************************************************
% Content: 
% In this file, the real-world Solgenia scenario is being evaluated using
% the elliptic cylinder, extruded superellipse and FCDS 3D EOT measurement
% models, as well as the SDFS, SH, and 3D GP comparison methods.
%
% The 3D GP model is taken from: https://github.com/Metu-Sensor-Fusion-Lab
%
% *************************************************************************

%% start of function code
close all
clear
clc

% random seed
rng(2,"twister")

% load data
load('data.mat')

% scenario measurements
meas = cellfun(@double,sol_meas,'UniformOutput',false);
nSamples = 300;

% simulation environment settings
do_plots = false;

 % choose if own thesis measurement models or comparison methods are
% evaluated. Choose either 'own', or 'comp'.
methods = 'own';

% set number measurements taken for a single update
n_upd = 20;

if strcmp(methods,'own')
    % filter settings. source either 'radial' or 'projected'. filter either 
    % 'ERHM' or 'GAM'.
    artificial_noise = true;
    filter = 'ERHM';
    source = 'radial';
elseif strcmp(methods,'comp')
    disp('Comparison methods are evaluated now!')
end

% set parameters for each measurement model applied in this simulation
set_parameters

%% start simulation
% start simulation
for k = 1:nSamples
    % display timestep
    disp(['timestep: ',num2str(k)])
    
    if strcmp(methods,'own')
        if strcmp(source,'projected')
            n_upd = size(meas{k,1},2);
            % elliptic cylinder UKF parameters
            nx = 9;
            if strcmp(filter,'ERHM')
                nxu = nx + n_upd;
            elseif strcmp(filter,'GAM')
                nxu = nx;
            end
            lambda_cyl = alpha_UKF^2*(nxu + kappa) - nxu;  
            % calculate weights for sigma points
            wm_cyl = zeros(1,2*nxu+1); wc_cyl = zeros(1,2*nxu+1);
            wm_cyl(1) = lambda_cyl/(nxu + lambda_cyl);
            wc_cyl(1) = lambda_cyl/(nxu + lambda_cyl) + (1 - alpha_UKF^2 + beta_UKF);
            wm_cyl(2:2*nxu + 1) = 1/(2*(nxu + lambda_cyl));
            wc_cyl(2:2*nxu + 1) = 1/(2*(nxu + lambda_cyl));
            % superellipse UKF parameters
            nx = 11;
            if strcmp(filter,'ERHM')
                nxu = nx + n_upd;
            elseif strcmp(filter,'GAM')
                nxu = nx;
            end
            lambda_super = alpha_UKF^2*(nxu + kappa) - nxu;
            % calculate weights for sigma points
            wm_super = zeros(1,2*nxu+1); wc_super = zeros(1,2*nxu+1);
            wm_super(1) = lambda_super/(nxu + lambda_super);
            wc_super(1) = lambda_super/(nxu + lambda_super) + (1 - alpha_UKF^2 + beta_UKF);
            wm_super(2:2*nxu + 1) = 1/(2*(nxu + lambda_super));
            wc_super(2:2*nxu + 1) = 1/(2*(nxu + lambda_super));
        end
    end

    if k == 1
        % do state initialization
        do_initialization
    else
        % do prediction steps 
        do_predictions
    end
    
    % generate random permutation of measurement set
    idxp = randperm(size(meas{k,1},2));

    % do measurement updates 
    do_updates

    % visualize results
    if do_plots
        plot_results
    end

    % calculate intersection over union using boat reference
    do_ious
end

%% evaluation
if strcmp(methods,'own')
    %% evaluate elliptic cylinder filter
    % constrained values
    X_cyl(7:9,:) = c1(X_cyl(7:9,:),0,'lower');

    % calculate rmses
    rmse_or_cyl = sqrt((X_cyl(5,:) - X_Ref(5,:)).^2);
    rmse_h_cyl = sqrt((X_cyl(9,:) - h_ref).^2);

    %% evaluate FCDS filter
    if strcmp(filter,'GAM')
        % constrained values
        X_fcds(7,:) = c1(X_fcds(7,:),0,'lower');
        % calculate rmses
        rmse_or_fcds = sqrt((X_fcds(5,:) - X_Ref(5,:)).^2);
        rmse_h_fcds = sqrt((X_fcds(7,:) - h_ref).^2);
    elseif strcmp(filter,'ERHM')
        % constrained values
        X_line(2,:) = c1(X_line(2,:),0,'lower');
        % calculate rmses
        rmse_or_fcds = sqrt((X_fcds(4,:) - X_Ref(5,:)).^2);
        rmse_h_fcds = sqrt((X_line(2,:) - h_ref).^2);
    end

    %% evaluate superellipse filter
    % constrained values
    X_super(7:9,:) = c1(X_super(7:9,:),0,'lower');
    X_super(10,:) = c1(X_super(10,:),1,'lower');
    X_super(11,:) = c2(X_super(11,:),-1,1);  

    % calculate rmses
    rmse_or_super = sqrt((X_super(5,:) - X_Ref(5,:)).^2);
    rmse_h_super = sqrt((X_super(9,:) - h_ref).^2);

    %% plot results
    % start figure
    tiledlayout(2,2)
    set(gcf,'Color','w')
    set(gcf,'Position',[2063,308,1074,456])
    
    % plot orientation rmse
    nexttile
    plot(1:nSamples,rmse_or_cyl,'LineWidth',2)
    hold on
    plot(1:nSamples,rmse_or_super,'LineWidth',2)
    plot(1:nSamples,rmse_or_fcds,'LineWidth',2)
    ylabel('RMSE','Interpreter','latex')
    title('\textbf{Orientation RMSE}','Interpreter','latex')
    % axes
    set(gca,'TickLength',[0 0])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    
    % plot height rmse
    nexttile
    plot(1:nSamples,rmse_h_cyl,'LineWidth',2)
    hold on
    plot(1:nSamples,rmse_h_super,'LineWidth',2)
    plot(1:nSamples,rmse_h_fcds,'LineWidth',2)
    ylabel('RMSE','Interpreter','latex')
    title('\textbf{Height RMSE}','Interpreter','latex')
    % axes
    set(gca,'TickLength',[0 0])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    
    % plot IoUs
    nexttile
    plot(1:nSamples,iou_cyl,'LineWidth',2)
    hold on
    plot(1:nSamples,iou_super,'LineWidth',2)
    plot(1:nSamples,iou_fcds,'LineWidth',2)
    plot(1:nSamples,iou_fcds_full,'LineWidth',2)
    xlabel('Time step','Interpreter','latex')
    ylabel('IoU','Interpreter','latex')
    title('\textbf{Intersection over Union}','Interpreter','latex')
    % axes
    set(gca,'TickLength',[0 0])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    
    % legend entries
    leg = legend('Elliptic Cylinder','Superellipse','FCDS','FCDS full');
    leg.Interpreter = "latex";
    leg.Orientation = "horizontal";
    leg.NumColumns = 5;
    leg.Box = "off";
    leg.Layout.Tile = "south";
    leg.FontSize = 12;
elseif strcmp(methods,'comp')
    %% evaluate SDFS filter
    % calculate rmses
    rmse_or_sdfs = sqrt((X_sdfs(5,:) - X_Ref(5,:)).^2);
    rmse_h_sdfs = sqrt((hs_sdfs - h_ref).^2);

    %% evaluate SH filter
    % calculate rmses
    rmse_or_sh = sqrt((X_sh(5,:) - X_Ref(5,:)).^2);
    rmse_h_sh = sqrt((hs_sh - h_ref).^2);

    %% evaluate GP filter
    % calculate orientation RMSE
    or_GP = zeros(1,nSamples);
    for i = 1:nSamples
        qq = quat2eul(quat_gp(:,i)');
        or_GP(i) = qq(1);
    end
    % calculate rmses
    rmse_or_gp = sqrt((or_GP - X_Ref(5,:)).^2);
    rmse_h_gp = sqrt((hs_gp - h_ref).^2);

    %% plot results
    % start figure
    tiledlayout(2,2)
    set(gcf,'Color','w')
    set(gcf,'Position',[2063,308,1074,456])
    
    % plot orientation rmse
    nexttile
    plot(1:nSamples,rmse_or_sdfs,'LineWidth',2)
    hold on
    plot(1:nSamples,rmse_or_sh,'LineWidth',2)
    plot(1:nSamples,rmse_or_gp,'LineWidth',2)
    ylabel('RMSE','Interpreter','latex')
    title('\textbf{Orientation RMSE}','Interpreter','latex')
    % axes
    set(gca,'TickLength',[0 0])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    
    % plot height rmse
    nexttile
    plot(1:nSamples,rmse_h_sdfs,'LineWidth',2)
    hold on
    plot(1:nSamples,rmse_h_sh,'LineWidth',2)
    plot(1:nSamples,rmse_h_gp,'LineWidth',2)
    ylabel('RMSE','Interpreter','latex')
    title('\textbf{Height RMSE}','Interpreter','latex')
    % axes
    set(gca,'TickLength',[0 0])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    
    % plot simplified IoUs
    nexttile
    plot(1:nSamples,iou_sdfs,'LineWidth',2)
    hold on
    plot(1:nSamples,iou_sh,'LineWidth',2)
    plot(1:nSamples,iou_gp,'LineWidth',2)
    xlabel('Time step','Interpreter','latex')
    ylabel('IoU','Interpreter','latex')
    title('\textbf{Simplified Intersection over Union}','Interpreter','latex')
    % axes
    set(gca,'TickLength',[0 0])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)

    % plot full IoUs
    nexttile
    plot(1:nSamples,iou_sdfs_full,'LineWidth',2)
    hold on
    plot(1:nSamples,iou_sh_full,'LineWidth',2)
    plot(1:nSamples,iou_gp_full,'LineWidth',2)
    xlabel('Time step','Interpreter','latex')
    ylabel('IoU','Interpreter','latex')
    title('\textbf{Full Intersection over Union}','Interpreter','latex')
    % axes
    set(gca,'TickLength',[0 0])
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    
    % legend entries
    leg = legend('SDFS','SH','GP');
    leg.Interpreter = "latex";
    leg.Orientation = "horizontal";
    leg.NumColumns = 5;
    leg.Box = "off";
    leg.Layout.Tile = "south";
    leg.FontSize = 12;
end

% end of function code
% *************************************************************************
%
%
% *************************************************************************
% end of document "main.m"
% *************************************************************************