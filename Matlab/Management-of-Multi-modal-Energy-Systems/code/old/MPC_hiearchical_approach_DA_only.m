%%

%
% This code runs the predictive multi-energy management system.

%
% This upper level MPC represents the continuous trading planning finding place
% shortly before the delivery and closing at 12:00. The numerical values are taken form the following paper
%
%   "Multi-stage optimal energy amangement of multi-energy microgrid in
%    deregulated electricity markets", Wang et al., Applied
%               Energy(2022).
%
% Copyright Institute of Automatic Control, Leibniz University Hannover.

%%
function [index_DA, index_CT, index_mon_DA, index_mon_CT] = MPC_hiearchical_approach_DA_only(T_DA,T_CT,days,noise, int_method)

clc
close all
yalmip('clear')


    %% load parameter  
    p = load_parameters();

    % Upper level MPC: Day ahead
    T_SC = 12;                           % scheduling time of the day ahead
    N_pred_DA = T_DA/p.Delta_DA;             
    N_SC = T_SC/p.Delta_DA;                 
    N_DA = N_pred_DA + N_SC;	          % discrete-time prediction horizon 
    N24 = 24;
    % Lower level MPC: Continuous trading
    N_CT = T_CT/p.Delta_CT;        % discrete-time prediction horizon [h]
    simStart = 1;
    sim_time = 24*days;               % simulate for x days in [h]
    simEnd = sim_time/p.Delta_CT;   % number of iterations

    %% Setup figures
    f1 = figure('name', 'Forecast uncertainty profiles for DA',...
                'Position',  [10, 100, 500, 500]);
    set(findall(gcf,'-property','FontSize'),'FontSize',20)

    f2 = figure('name', '24h-1h Day ahead simulations (next day)',...
                'Position',  [110, 100, 1000, 800]);
    set(findall(gcf,'-property','FontSize'),'FontSize',20)

    f3 = figure('name', '24h-1h Closed-loop simulations (next day)',...
                'Position',  [210, 100, 1000, 800]);
    set(findall(gcf,'-property','FontSize'),'FontSize',20)


   %% load measurement data
    
   % load real data
   load ('../data/8760h-01p/DA_real.mat');
   load ('../data/8760h-01p/CT_real.mat');

    [out_DA, out_CT] = load_data(T_DA, T_CT, noise);  % no noise
  
    n_DA = length(out_DA);
    n_CT = length(out_CT);
    %out_DA = [PL_DA; PV_DA; WT_DA;HL_DA;v_DA; u_DA];   
    %out_CT = [PL_CT; PV_CT; WT_CT;HL_CT;v_CT; u_CT];   

     if noise ==0

     PL_DA = out_DA(1,:);
     PV_DA = out_DA(2,:);
     WT_DA = out_DA(3,:);
     HL_DA = out_DA(4,:);
     v_DA = out_DA(5,:);
     u_DA = out_DA(6,:);
    % 
      
     PL_CT = out_CT(1,:);
     PV_CT = out_CT(2,:);
     WT_CT = out_CT(3,:);
     HL_CT = out_CT(4,:);
     v_CT = out_CT(5,:);
     u_CT = out_CT(6,:);

     else
        
     PL_DA = out_DA(1:n_DA/6,:);
     PV_DA = out_DA(n_DA/6+1:2*n_DA/6,:);
     WT_DA = out_DA(2*n_DA/6+1:3*n_DA/6,:);
     HL_DA = out_DA(3*n_DA/6+1:4*n_DA/6,:);
     v_DA = out_DA(4*n_DA/6+1:5*n_DA/6,:);
     u_DA = out_DA(5*n_DA/6+1:end,:);
    % 
      
     PL_CT = out_CT(1:n_CT/6,:);
     PV_CT = out_CT(n_CT/6+1:2*n_CT/6,:);
     WT_CT = out_CT(2*n_CT/6+1:3*n_CT/6,:);
     HL_CT = out_CT(3*n_CT/6+1:4*n_CT/6,:);
     v_CT = out_CT(4*n_CT/6+1:5*n_CT/6,:);
     u_CT = out_CT(5*n_CT/6+1:end,:);
    
     end

  
   
    %% Plot load and generation profiles
    % uncertainties_DA = [PL_DA(times); HL_DA(times); PV_DA(times); WT_DA(times); u_DA(times); v_DA(times)];
    set (0, 'CurrentFigure', f1);
    
    ax31 = subplot(4,1,1);
    plot(v_DA(1:simEnd/4+N_DA))
    hold on
    plot(u_DA(1:simEnd/4+N_DA))
    grid on
    title('Prices')
    xlabel('Step in (every 1h)')
    ylabel('Price in dollars')
    xlim([1,simEnd/4+N_DA])

    ax32 = subplot(4,1,2);
    plot(PL_DA(1:simEnd/4+N_DA));
    grid on
    hold on
    plot(HL_DA(1:simEnd/4+N_DA));
    title('Load profile')
    xlabel('Step in (every 1h)')
    ylabel('Power in kW')
    xlim([1,simEnd/4+N_DA])

    ax33 = subplot(4,1,3);
    plot(PV_DA(1:simEnd/4+N_DA));
    grid on
    title('PV generation')
    xlabel('Step in (every 1h)')
    ylabel('Power in kW')
    xlim([1,simEnd/4+N_DA])

    ax34 = subplot(4,1,4);
    plot(WT_DA(1:simEnd/4+N_DA));
    grid on
    title('Wind Turbine generation')
    xlabel('Step in (every 1h)')
    ylabel('Power in kW')
    xlim([1,simEnd/4+N_DA])

    drawnow;



%% MPC Closed loop simulations
    [var_DA, state_DA] = generate_vars(N_DA);
    [var_CT, state_CT] = generate_vars(N_CT);
    % DA Initial conditions measured 12h ago
    x_DA_measure = [0.3;600];
    
    x = [];
    u = [];
    Pg_DA_Cls = [];
    x_OL_vec = [];
    u_OL_vec = [];
    Pg_DA_vec = [];
    P_CHP_DA_vec = [];
    Pgrid_Cls = [];
    P_CHP_Cls = [];
   % H_GB_Cls = [];
    index_DA = 0;
    index_CT = 0;
    index_mon_DA = 0;
    index_mon_CT = 0;
    zz = 49;
    rr = 1;
    s = 0;
   % Receding horizon control
    tic
for ii = simStart:simEnd
    clc
    %% Step I: Day ahead 
    % assume we are at 12h one day before delivery, every 24h iterations, the upper level is triggered
   if ii == rr 
    % update uncertainty profiles from the forecast for the next N=24h hours
    kk = 0;
    s = s+1;  % number of times DA is getriggered
    % take the measurement of (12:00) of the previous day
    
    SoC_init_DA = x_DA_measure(1);
    H_init_DA = x_DA_measure(2);

    times_DA = ceil(ii/4):N_DA-1+ceil(ii/4);  % take hourly measurements
    if noise ==0
    uncertainties_DA = [PL_DA(times_DA); HL_DA(times_DA); PV_DA(times_DA); ...
        WT_DA(times_DA); u_DA(times_DA); v_DA(times_DA)];
    else
     uncertainties_DA = [PL_DA(s,:); HL_DA(s,:); PV_DA(s,:); ...
       WT_DA(s,:); u_DA(s,:); v_DA(s,:)];
    end
    % setup mpc for day ahead
    mpc_DA = setup_day_ahead(var_DA, state_DA, uncertainties_DA, N_DA, p);
    % re-Initialisation
    init_DA = {SoC_init_DA, H_init_DA};
    [sol_DA,flag_DA] = mpc_DA(init_DA);
    % out = {SoC,H,Pc,Pd,Hc,Hd,P_EB,P_CHP,I_ESS,I_HSS,H_EB,H_GB,Pgrid,cost_DA};
    if flag_DA > 0
        fprintf('Exitflag DA is %3d',flag_DA)
       break
    end
    % states starting from the following day at midnight
    next_day = N_SC+1:N24+N_SC;
    x_OL(1,:) = sol_DA{1}(N_SC+1:N24+N_SC+1); % SOC
    x_OL(2,:) = sol_DA{2}(N_SC+1:N24+N_SC+1); % H
    % inputs
    u_OL(1,:) = sol_DA{3}(next_day); % Pc
    u_OL(2,:) = sol_DA{4}(next_day); % Pd
    u_OL(3,:) = sol_DA{5}(next_day); % Hc
    u_OL(4,:) = sol_DA{6}(next_day); % Hd
    Pg_DA = sol_DA{13}(next_day); % Pg
    P_CHP_DA = sol_DA{8}(next_day); % PCHP

    input_DA = u_OL(1:4,:); % input to be used by the MES later

    P_DA_var = [sol_DA{3}(next_day); sol_DA{4}(next_day); Pg_DA; sol_DA{8}(next_day)];
    H_DA_var = [sol_DA{5}(next_day); sol_DA{6}(next_day); sol_DA{12}(next_day)];
    mon_var = [sol_DA{13}(next_day); sol_DA{8}(next_day); sol_DA{12}(next_day)];

    % DA performance index
    times_P = ceil(ii/4)+N_SC:N_pred_DA-1+ceil(ii/4)+N_SC;  % start on the next day (i.e., after 12h)
    index_DA = index_DA + setup_cost_DA(v_DA_real(times_P),u_DA_real(times_P),P_DA_var,H_DA_var, N24, p);
    index_mon_DA = index_mon_DA + ...
    monetary_index_DA(v_DA_real(times_P),u_DA_real(times_P),mon_var, p);


    %% interpolate DA solution 
   
     Pg_DA_rshp = [];
     input_DA_rshp = [];
            for jj = 1:size(Pg_DA,2)
                Pg_DA_rshp = [Pg_DA_rshp repmat(Pg_DA(jj),1,p.Delta_DA/p.Delta_CT)];
                input_DA_rshp = [input_DA_rshp repmat(input_DA(:,jj),1,p.Delta_DA/p.Delta_CT)];
            end
    Pg_DA_rshp = [Pg_DA_rshp Pg_DA_rshp(end)*ones(1,p.Delta_DA*T_CT/p.Delta_CT)]; % for the horizon N_CT
    
    if (int_method == "m2" && T_DA == 48)
          Pg_DA = sol_DA{13}(N_SC+1:end); % only
          Pg_DA_rshp = [];
            for jj = 1:size(Pg_DA,2)
                Pg_DA_rshp = [Pg_DA_rshp repmat(Pg_DA(jj),1,p.Delta_DA/p.Delta_CT)];
            end
    else 
        if (int_method == "m2" && T_DA ~= 48)
        fprintf("ERROR: DA horizon is not 48 \n")
        break;
        end
    end


        rr = rr + p.Delta_DA*N24/p.Delta_CT;  % trigger DA every 24h 
       
        x_OL_vec = [x_OL_vec, x_OL];
        u_OL_vec = [u_OL_vec, u_OL];
        Pg_DA_vec = [Pg_DA_vec, Pg_DA(1:N24)];
        P_CHP_DA_vec = [P_CHP_DA_vec, P_CHP_DA];

        

      % plots of the Day ahead
        set (0, 'CurrentFigure', f2);

        subplot(2,4,1), cla
        hold on
        plot(1:N24*s+s, x_OL_vec(1,:),'b')
        grid on
        ylabel('SoC (in p.u.)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])
        drawnow

        subplot(2,4,2), cla
        hold on
        plot(1:N24*s+s, x_OL_vec(2,:),'r')
        grid on
        ylabel('H (in kWh)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])
        drawnow

        subplot(2,4,3), cla
        hold on
        plot(1:N24*s, u_OL_vec(1,:), 'b')
        grid on
        ylabel('Pc (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])
        drawnow

        subplot(2,4,4), cla
        hold on
        plot(1:N24*s, u_OL_vec(2,:), 'b')
        grid on
        ylabel('Pd (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])
        drawnow

        subplot(2,4,5), cla
        hold on
        plot(1:N24*s, u_OL_vec(3,:),'r')
        grid on
        ylabel('Hc (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])
        drawnow

        subplot(2,4,6), cla
        hold on
        plot(1:N24*s, u_OL_vec(4,:),'r')
        grid on
        ylabel('Hd (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])
        drawnow

        subplot(2,4,7), cla
        hold on
        plot(1:N24*s, Pg_DA_vec, 'g')
        hold on
        grid on
        ylabel('Pgrid (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])
        drawnow

        subplot(2,4,8), cla
        hold on
        plot(1:N24*s, P_CHP_DA_vec, 'g')
        hold on
        grid on
        ylabel('P_{CHP} (in kW)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])
        drawnow


   end



     %% Step II: Continuous trading (shortly before midnight)
        if ii ==1
            x_CT_measure = x_OL(:,1);  % take the predicted state by the DA
        end
        SoC_init_CT = x_CT_measure(1);
        H_init_CT = x_CT_measure(2);

        % update uncertainty profiles from the CT forecast
        times = ii:N_CT-1+ii;
        kk = kk+1;
        if noise ==0
        uncertainties_CT = [PL_CT(times); HL_CT(times); PV_CT(times); ...
        WT_CT(times); u_CT(times); v_CT(times)];
        else
        uncertainties_CT = [PL_CT(ii,:); HL_CT(ii,:); PV_CT(ii,:); ...
        WT_CT(ii,:); u_CT(ii,:); v_CT(ii,:)];
        end

        % STEP 3: solve CT mpc problem
        mpc_CT = setup_continuous_trading(var_CT,state_CT, uncertainties_CT, Pg_DA_rshp(kk:N_CT-1+kk),N_CT, p);
        % re-Initialisation
        init_CT = {SoC_init_CT, H_init_CT};
        [sol_CT,flag_CT] = mpc_CT(init_CT);
        if flag_CT > 0
            fprintf('Exitflag CT is %3d',flag_CT)
            break;
        end

        % STEP 4: update the measurements
        % out = {SoC,H,Pc,Pd,Hc,Hd,P_EB,P_CHP,I_ESS,I_HSS,H_EB,H_GB(12),Pgrid(13),cost_CT};
        input_CT = [sol_CT{3}(1);sol_CT{4}(1); sol_CT{5}(1); sol_CT{6}(1)];  % only first input in the sequence
        P_CT_var = [sol_CT{3}(1); sol_CT{4}(1); sol_CT{13}(1); sol_CT{8}(1)];
        H_CT_var = [sol_CT{5}(1); sol_CT{6}(1); sol_CT{12}(1)];
        vec = [sol_CT{13}(1); sol_CT{8}(1); sol_CT{12}(1)];
        %t = [t, Delta_CT];
        Pgrid_Cls = [Pgrid_Cls, sol_CT{13}(1)];
        P_CHP_Cls = [P_CHP_Cls, sol_CT{8}(1)];
        % H_GB_Cls = [H_GB_Cls, sol_CT{12}(1)];
        Pg_DA_Cls = [Pg_DA_Cls, Pg_DA_rshp(kk)];

        % CT performance index setup_cost_CT(v_GAS_CT,v,u,P_var,H_var,Pg_DA_rshp, N_CT, Delta_CT)
        index_CT = index_CT + setup_cost_CT(v_CT_real(ii),u_CT_real(ii),P_CT_var,H_CT_var, Pg_DA_rshp(kk), 1, p);
        index_mon_CT = index_mon_CT + ...
        monetary_index_CT(v_CT_real(ii),u_CT_real(ii),vec, Pg_DA_rshp(kk), p);

 
        [out1, out2] = dynamics(x_CT_measure, input_CT, p); % simulate the closed-loop system
        x_CT_measure = [out1; out2];
         
        x = [x, x_CT_measure];
        u = [u, input_CT];
        
        if ii == zz  % update DA measurements at 12h
            x_DA_measure = x_CT_measure;
            zz = zz + p.Delta_DA*T_SC/p.Delta_CT; % save the x measurement every 12h
        end
       

        % plots
        set (0, 'CurrentFigure', f3);

        subplot(2,4,1), cla
        hold on
        plot(x(1,:),'bo-')
        grid on
        ylabel('SoC (in p.u.)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,2), cla
        hold on
        plot(x(2,:),'ro-')
        grid on
        ylabel('H (in kWh)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,3), cla
        hold on
        plot(u(1,:),'bo-')
        grid on
        ylabel('Pc (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,4), cla
        hold on
        plot(u(2,:),'bo-')
        grid on
        ylabel('Pd (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,5), cla
        hold on
        plot(u(3,:),'ro-')
        grid on
        ylabel('Hc (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,6), cla
        hold on
        plot(u(4,:),'ro-')
        grid on
        ylabel('Hd (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,7), cla
        hold on
        plot(Pgrid_Cls,'bo-')
        hold on 
        plot(Pg_DA_Cls,'go-')
        grid on
        ylabel('Pgrid_{CT}/Pgrid_{DA} (in kW)' )
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,4,8), cla
        hold on
        plot(P_CHP_Cls,'go-')
        grid on
        ylabel('P_{CHP} (in kW)')
        xlabel('Steps (in every 15 min)')
        drawnow

       

    
 end
 toc


end
