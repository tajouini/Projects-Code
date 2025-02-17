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
function [index_DA, index_CT, index_mon_DA, index_mon_CT, flag_DA, flag_CT] = MPC_hiearchical_approach(T_DA,T_CT,days,noise)

clc
close all
yalmip('clear')
%addpath('figs/Matlab2tikz')

    %% load parameter  
    p = load_parameters();

    % Upper level MPC: Day ahead
    T_SC = 12;                           % scheduling time of the day ahead
    N_pred_DA = T_DA/p.Delta_DA;             
    N_SC = T_SC/p.Delta_DA;                 
    N_DA = N_pred_DA + N_SC;	          % discrete-time prediction horizon 
    N24 = 24;
    % Lower level MPC: Continuous trading
    T_0 = 1/4;
    N_pred_CT = T_CT/p.Delta_CT;        
    N_CT = N_pred_CT + T_0/p.Delta_CT; % discrete-time prediction horizon [h]
    simStart = 1;
    sim_time = 24*days;               % simulate for x days in [h]
    simEnd = sim_time/p.Delta_CT;   % number of iterations

    
    %% Setup figures
    f1 = figure('name', 'Forecast uncertainty profiles for DA',...
                'Position', [100, -90, 500, 500]);
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    f2 = figure('name', 'Day ahead simulations (next day)',...
                'Position',  [200, -90, 1000, 800]);
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
    f3 = figure('name', 'Closed-loop simulations (next day)',...
                'Position',   [100, -90, 1000, 800]);
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
   % f4 = figure('name', 'Open vs. closed-loop simulations (next day)',...
           %     'Position',   [100, -90, 1000, 800]);
   % set(findall(gcf,'-property','FontSize'),'FontSize',20)


   %% load measurement data
    
   % load real data
   load ('../data/8760h-01p/DA_real.mat');
   load ('../data/8760h-01p/CT_real.mat');
   load ('../data/8760h-01p/noise.mat');

    [out_DA, out_CT] = load_data(T_DA, T_CT, noise);  % no noise
  
    n_DA = length(out_DA);
    n_CT = length(out_CT);
    % out_DA = [PL_DA; PV_DA; WT_DA;HL_DA;v_DA; u_DA];   
    % out_CT = [PL_CT; PV_CT; WT_CT;HL_CT;v_CT; u_CT];   

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
    hold on 
    plot(WT_DA(1:simEnd/4+N_DA));
    grid on
    title('PV and WT generation')
    xlabel('Step in (every 1h)')
    ylabel('Power in kW')
    xlim([1,simEnd/4+N_DA])

    % ax34 = subplot(4,1,4);
    % plot(WT_DA(1:simEnd/4+N_DA));
    % grid on
    % title('Wind Turbine generation')
    % xlabel('Step in (every 1h)')
    % ylabel('Power in kW')
    % xlim([1,simEnd/4+N_DA])

    % cleanfigure;
    % matlab2tikz('real-data.tex', 'parseStrings',false)


     
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

  
   
  %% MPC Closed loop simulations
    [var_DA, state_DA] = generate_vars(N_DA);
    [var_CT, state_CT] = generate_vars(N_CT);
    % DA Initial conditions measured 12h ago
    x_DA_measure = p.x0;

    % day ahead inputs and states
    x_OL_vec = [];
    u_OL_vec = [];
    Pg_DA_vec = [];
    P_CHP_DA_vec = [];
    deltaP_DA_vec = [];

    % closed-loop inputs and states
    x = [];
    u = [];
    Pg_DA_Cls = [];
    x_DA_Cls = [];
    Pgrid_Cls = [];
    P_CHP_Cls = [];
    H_GB_Cls = [];

    % CT interpolated day-ahead inputs and states 
    Pc_DA_int_Cls = [];
    Pd_DA_int_Cls = [];
    Hc_DA_int_Cls = [];
    Hd_DA_int_Cls = [];
    H_GB_DA_int_Cls = [];
    P_CHP_DA_int_Cls = [];
    deltaP_CT_vec = [];

   
    % initialize indices
    index_DA = 0;
    index_CT = 0;
    index_mon_DA = 0;
    index_mon_CT = 0;
    zz = 49;
    rr = 1;
    s = 0;
    %flag_CT = 0;
    %factor = 0.02;  % level of noise corruption
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
    % states starting from one hour before midnight
    next_day = N_SC:N_pred_DA+N_SC; % forward the whole sequence of length T_DA
    % for CT interpolation we also take one hour before midnight
    x_OL(1,:) = sol_DA{1}(N_SC:N_pred_DA+N_SC+1); % SOC
    x_OL(2,:) = sol_DA{2}(N_SC:N_pred_DA+N_SC+1); % H

   
    % inputs
    u_OL(1,:) = sol_DA{3}(next_day); % Pc
    u_OL(2,:) = sol_DA{4}(next_day); % Pd
    u_OL(3,:) = sol_DA{5}(next_day); % Hc
    u_OL(4,:) = sol_DA{6}(next_day); % Hd


    Pg_DA = sol_DA{13}(next_day); % Pg
    P_CHP_DA = sol_DA{8}(next_day); % PCHP
    H_EB_DA = sol_DA{11}(next_day); % H_EB
    H_GB_DA = sol_DA{12}(next_day); % H_GB

  
    % DA performance index 
    if ii>1   % discared the first day
    
    next_day_only = N_SC+1:N_pred_DA+N_SC;
    P_DA_var = [sol_DA{3}(next_day_only); sol_DA{4}(next_day_only); sol_DA{13}(next_day_only); sol_DA{8}(next_day_only)];
    H_DA_var = [sol_DA{5}(next_day_only); sol_DA{6}(next_day_only); sol_DA{12}(next_day_only)];
    mon_var = [sol_DA{13}(next_day_only); sol_DA{8}(next_day_only); sol_DA{12}(next_day_only)];

    times_P = ceil(ii/4)+N_SC:N24-1+ceil(ii/4)+N_SC;  % start on the next day (i.e., after 12h)
    index_DA = index_DA + setup_cost_DA(v_DA_real(times_P),u_DA_real(times_P),P_DA_var,H_DA_var, N24, p);
    index_mon_DA = index_mon_DA + ...
    monetary_index_DA(v_DA_real(times_P),u_DA_real(times_P),mon_var, p);

    end

    %% interpolate DA solution 
    
     Pg_DA_rshp = [];
     x_DA_int = [];
     Pc_DA_int = [];
     Pd_DA_int = [];
     Hc_DA_int = [];
     Hd_DA_int = [];
     P_CHP_DA_int = [];
     H_EB_DA_int = [];
     H_GB_DA_int = [];
     
            for jj = 1:size(Pg_DA,2)
                Pg_DA_rshp = [Pg_DA_rshp repmat(Pg_DA(jj),1,p.Delta_DA/p.Delta_CT)];
                x_DA_int = [x_DA_int  interpolate_pts(x_OL(:,jj),u_OL(:,jj),p)];
                % interpolate inputs
                Pc_DA_int = [Pc_DA_int repmat(u_OL(1,jj),1,p.Delta_DA/p.Delta_CT)];
                Pd_DA_int = [Pd_DA_int repmat(u_OL(2,jj),1,p.Delta_DA/p.Delta_CT)];
                Hc_DA_int = [Hc_DA_int repmat(u_OL(3,jj),1,p.Delta_DA/p.Delta_CT)];
                Hd_DA_int = [Hd_DA_int repmat(u_OL(4,jj),1,p.Delta_DA/p.Delta_CT)];
                P_CHP_DA_int = [P_CHP_DA_int repmat(P_CHP_DA(jj),1,p.Delta_DA/p.Delta_CT)];
                H_EB_DA_int = [H_EB_DA_int repmat(H_EB_DA(jj),1,p.Delta_DA/p.Delta_CT)];
                H_GB_DA_int = [H_GB_DA_int repmat(H_GB_DA(jj),1,p.Delta_DA/p.Delta_CT)];

            end
        Pg_DA_rshp = [Pg_DA_rshp Pg_DA_rshp(end)*ones(1,p.Delta_DA*T_CT/p.Delta_CT)]; % for the horizon N_CT
  

        % Remove the first 3 values because CT starts 15 min before midnight
        Pg_DA_rshp = Pg_DA_rshp(4:end);
        Pc_DA_int = Pc_DA_int(4:end);
        Pd_DA_int = Pd_DA_int(4:end);
        Hc_DA_int = Hc_DA_int(4:end);
        Hd_DA_int = Hd_DA_int(4:end);
        P_CHP_DA_int = P_CHP_DA_int(4:end);
        H_GB_DA_int = H_GB_DA_int(4:end);
        H_EB_DA_int = H_EB_DA_int(4:end);
        
        
        x_DA_int = x_DA_int(:,4:end); % remove the first three values

        rr = rr + p.Delta_DA*N24/p.Delta_CT;  % trigger DA every 24h 
     
        x_OL_vec = [x_OL_vec, x_OL(:,1:N24+2)]; % 25+1 
        u_OL_vec = [u_OL_vec, u_OL(:,1:N24+1)];
        Pg_DA_vec = [Pg_DA_vec, Pg_DA(1:N24+1)];
        P_CHP_DA_vec = [P_CHP_DA_vec, P_CHP_DA(1:N24+1)];
        
        deltaP_DA = uncertainties_DA(3,N_SC:end) + uncertainties_DA(4,N_SC:end)-...
            uncertainties_DA(1,N_SC:end);
        deltaP_DA_vec = [deltaP_DA_vec deltaP_DA];

      % plots of the Day ahead
     
        set (0, 'CurrentFigure', f2);

        subplot(2,5,1), cla
        hold on
        plot(x_OL_vec(1,:),'b')
        grid on
        ylabel('SoC (in p.u.)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])
        drawnow
        hold on

        set (0, 'CurrentFigure', f2);
        subplot(2,5,2), cla
        hold on
        plot(x_OL_vec(2,:),'r')
        grid on
        ylabel('H (in kWh)')
        xlabel('Steps (in every 1h)')
        xlim([1, simEnd*p.Delta_CT])
        drawnow

        
        set (0, 'CurrentFigure', f2);
         subplot(2,5,3), cla
         hold on
         plot(u_OL_vec(1,:), 'b')
         grid on
         ylabel('Pc (in kW)')
         xlabel('Steps (in every 1h)')
         xlim([1, simEnd*p.Delta_CT])
         drawnow

         set (0, 'CurrentFigure', f2);
         subplot(2,5,4), cla
         hold on
         plot(u_OL_vec(2,:), 'b')
         grid on
         ylabel('Pd (in kW)')
         xlabel('Steps (in every 1h)')
         xlim([1, simEnd*p.Delta_CT])
         drawnow

         set (0, 'CurrentFigure', f2);
         subplot(2,5,5), cla
         hold on
         plot(u_OL_vec(3,:),'r')
         grid on
         ylabel('Hc (in kW)')
         xlabel('Steps (in every 1h)')
         xlim([1, simEnd*p.Delta_CT])
         drawnow

         set (0, 'CurrentFigure', f2);
         subplot(2,5,6), cla
         hold on
         plot(u_OL_vec(4,:),'r')
         grid on
         ylabel('Hd (in kW)')
         xlabel('Steps (in every 1h)')
         xlim([1, simEnd*p.Delta_CT])
         drawnow

         set (0, 'CurrentFigure', f2);
         subplot(2,5,7), cla
         hold on
         plot(Pg_DA_vec, 'g')
         hold on
         grid on
         ylabel('Pgrid (in kW)')
         xlabel('Steps (in every 1h)')
         xlim([1, simEnd*p.Delta_CT])
         drawnow

         set (0, 'CurrentFigure', f2);
         subplot(2,5,8), cla
         hold on
         plot(P_CHP_DA_vec, 'g')
         hold on
         grid on
         ylabel('P_{CHP} (in kW)')
         xlabel('Steps (in every 1h)')
         xlim([1, simEnd*p.Delta_CT])
         drawnow

         set (0, 'CurrentFigure', f2);
         subplot(2,5,9), cla
         hold on
         plot(deltaP_DA_vec, 'g')
         hold on
         grid on
         ylabel('P_{WT}+P_{PV}-P_L (in kW)')
         xlabel('Steps (in every 1h)')
         xlim([1, simEnd*p.Delta_CT])
         drawnow

       


   end



     %% Step II: Continuous trading start at 15 min before delivery

        if ii ==1 % 23:45 measurement

          x_CT_prev = x_DA_int(:,1);
     
        end
     
        SoC_init_CT = x_CT_prev(1);
        H_init_CT = x_CT_prev(2);

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
        % For ii= 1 find an IC of the multi-energy system
        if ii ==1
        %mpc_CT = setup_continuous_trading(var_CT,state_CT, uncertainties_CT, Pg_DA_rshp(kk:N_CT-1+kk),N_CT, p);
        % we fix the input to the current one and optimize starting from the second input 
        sol_DA_int{3} = Pc_DA_int;
        sol_DA_int{4} = Pd_DA_int;
        sol_DA_int{5} = Hc_DA_int;
        sol_DA_int{6} = Hd_DA_int;
        sol_CT_prev = sol_DA_int;
        
        mpc_CT1 = setup_continuous_trading_fix_u1(sol_CT_prev,var_CT,state_CT, uncertainties_CT, Pg_DA_rshp(kk:N_CT-1+kk),N_CT, p);
        % re-Initialisation
        init_CT = {SoC_init_CT, H_init_CT};
        [sol_CT,flag_CT] = mpc_CT1(init_CT);
       

        if flag_CT > 0
            fprintf('Exitflag CT is %3d',flag_CT)
            break;
        end
        
        
       
        % take the second predicted state by CT as IC of the mmes to ensure feasibility
         x_CT_measure = [sol_CT{1}(2); sol_CT{2}(2)];%+0.01*[sol_CT{1}(2); sol_CT{2}(2)].*noise1(1:2)';
        
        else

        sol_CT_prev = sol_CT;
        mpc_CT_fix_u1 = setup_continuous_trading_fix_u1(sol_CT_prev, var_CT,state_CT, uncertainties_CT, Pg_DA_rshp(kk:N_CT-1+kk),N_CT, p);
        % re-Initialisation
        init_CT = {SoC_init_CT, H_init_CT};
        [sol_CT,flag_CT] = mpc_CT_fix_u1(init_CT);
       
        if flag_CT > 0
            fprintf('Exitflag CT is %3d',flag_CT)
            break;
        end
   
        end


        % take the second input in the sequence
        input_CT = [sol_CT{3}(2);sol_CT{4}(2); sol_CT{5}(2); sol_CT{6}(2)]; 

        % STEP 4: update the measurements
        P_CT_var = [sol_CT{3}(2); sol_CT{4}(2); sol_CT{13}(2); sol_CT{8}(2)];
        %P_CT_var = [sol_CT{3}(2:N_CT); sol_CT{4}(2:N_CT); sol_CT{13}(2:N_CT); sol_CT{8}(2:N_CT)];
        H_CT_var = [sol_CT{5}(2); sol_CT{6}(2); sol_CT{12}(2)];
        %H_CT_var = [sol_CT{5}(2:N_CT); sol_CT{6}(2:N_CT); sol_CT{12}(2:N_CT)];
        vec = [sol_CT{13}(2); sol_CT{8}(2); sol_CT{12}(2)];
        %vec = [sol_CT{13}(2:N_CT); sol_CT{8}(2:N_CT); sol_CT{12}(2:N_CT)];
        %t = [t, Delta_CT];
        Pgrid_Cls = [Pgrid_Cls, sol_CT{13}(2)];
        P_CHP_Cls = [P_CHP_Cls, sol_CT{8}(2)];
        H_GB_Cls = [H_GB_Cls, sol_CT{12}(2)];
        Pg_DA_Cls = [Pg_DA_Cls, Pg_DA_rshp(kk+1)];

        % the predicted SOC and heat by DA: verify indices
        x_DA_Cls = [x_DA_Cls, x_DA_int(:,kk+1)];  
        Pc_DA_int_Cls = [Pc_DA_int_Cls Pc_DA_int(kk+1)];
        Pd_DA_int_Cls = [Pd_DA_int_Cls Pd_DA_int(kk+1)];
        Hc_DA_int_Cls = [Hc_DA_int_Cls Hc_DA_int(kk+1)];
        Hd_DA_int_Cls = [Hd_DA_int_Cls Hd_DA_int(kk+1)];
        P_CHP_DA_int_Cls = [P_CHP_DA_int_Cls P_CHP_DA_int(kk+1)];
        H_GB_DA_int_Cls = [H_GB_DA_int_Cls H_GB_DA_int(kk+1)];
        
        deltaP_CT = uncertainties_CT(3,:) + uncertainties_CT(4,:)-...
        uncertainties_CT(1,:);
        deltaP_CT_vec = [deltaP_CT_vec deltaP_CT];

        % CT performance index setup_cost_CT(v_GAS_CT,v,u,P_var,H_var,Pg_DA_rshp, N_CT, Delta_CT)
        % the costs of the next 15 min (kk+1)
       if ii>94  % discared first day 
         index_CT = index_CT + setup_cost_CT(v_CT_real(ii+1),u_CT_real(ii+1),P_CT_var,H_CT_var, Pg_DA_rshp(kk+1), 1, p);
         index_mon_CT = index_mon_CT + ...
         monetary_index_CT(v_CT_real(ii+1),u_CT_real(ii+1),vec, Pg_DA_rshp(kk+1), 1, p);
       end

        % index_CT = index_CT + setup_cost_CT(v_CT_real(2:N_CT),u_CT_real(2:N_CT),P_CT_var,H_CT_var, Pg_DA_rshp(2:N_CT), N_CT-1, p);
        % index_mon_CT = index_mon_CT + ...
        % monetary_index_CT(v_CT_real(2:N_CT),u_CT_real(2:N_CT),vec, Pg_DA_rshp(2:N_CT), N_CT-1, p);
        % 

        % simulate the closed-loop system
        [out1, out2] = dynamics(x_CT_measure, input_CT, p); 
       
        % update the measurements
        x_CT_prev = x_CT_measure; % store 15 min old measurement
        x_CT_measure = [out1; out2]; % update the measurement
        x = [x, x_CT_measure];
        u = [u, input_CT];
        
        % update DA measurements every 24h
        if ii == zz  
        x_DA_measure =  x_CT_measure;
        zz = zz + p.Delta_DA*N24/p.Delta_CT; % save the x DA measurement at 12h
        end
       

        %% plots of closed-loop simulations
       
        set (0, 'CurrentFigure', f3);
        % 
        subplot(2,5,1), cla
        plot(x(1,:),'bo-')
        hold on
        plot(x_DA_Cls(1,:),'go-')
        grid on
        ylabel('SoC (in p.u.)')
        xlabel('Steps (in every 15 min)')
        drawnow

        subplot(2,5,2), cla
        plot(x(2,:),'bo-')
        hold on 
        plot(x_DA_Cls(2,:),'go-')
        grid on
        ylabel('H (in kWh)')
        xlabel('Steps (in every 15 min)')
        drawnow




         subplot(2,5,3), cla
         plot(u(1,:),'bo-')
         hold on 
         plot(Pc_DA_int_Cls,'go-')
         grid on
         ylabel('Pc (in kW)')
         xlabel('Steps (in every 15 min)')
         drawnow
         % 
         subplot(2,5,4), cla
         plot(u(2,:),'bo-')
         hold on 
         plot(Pd_DA_int_Cls,'go-')
         grid on
         ylabel('Pd (in kW)')
         xlabel('Steps (in every 15 min)')
         drawnow

         % 
         subplot(2,5,5), cla
         plot(u(3,:),'bo-')
         hold on
         plot(Hc_DA_int_Cls,'go-')
         grid on
         ylabel('Hc (in kW)')
         xlabel('Steps (in every 15 min)')
         drawnow

         subplot(2,5,6), cla
         plot(u(4,:),'bo-')
         hold on
         plot(Hd_DA_int_Cls,'go-')
         grid on
         ylabel('Hd (in kW)')
         xlabel('Steps (in every 15 min)')
         drawnow
         % 
         subplot(2,5,7), cla
         hold on
         plot(Pgrid_Cls,'b-')
         hold on 
         plot(Pg_DA_Cls,'g-')
         grid on
         ylabel('Pgrid_{CT}/Pgrid_{DA} (in kW)' )
         xlabel('Steps (in every 15 min)')
         drawnow
         % 
         subplot(2,5,8), cla
         plot(P_CHP_Cls,'bo-')
         hold on
         plot(P_CHP_DA_int_Cls,'go-')
         grid on
         ylabel('P_{CHP} (in kW)')
         xlabel('Steps (in every 15 min)')
         drawnow

         % if ii == simEnd
         %     cleanfigure;
         %     matlab2tikz('Pgrid4824.tex', 'parseStrings',false)
         % end

         % subplot(2,5,9), cla
         % plot(sol_CT{11},'bo-')
         % hold on
         % plot(H_EB_Cls,'go-')
         % plot(H_EB_DA_int,'go-')
         % grid on
         % ylabel('H_{EB} (in kW)')
         % xlabel('Steps (in every 15 min)')
         % drawnow
         % 
         % subplot(2,5,9), cla
         % plot(H_GB_Cls,'bo-')
         % hold on 
         % plot(H_GB_DA_int_Cls,'go-')
         % grid on
         % ylabel('H_{GB} (in kW)')
         % xlabel('Steps (in every 15 min)')
         % drawnow

        

        %% open-loop simulations
        % set (0, 'CurrentFigure', f4);
        % 
        % subplot(2,6,1), cla
        % plot(sol_CT{1}(:),'bo-')
        % hold on
        % plot(x_DA_int(1,:),'go-')
        % grid on
        % ylabel('SoC (in p.u.)')
        % xlabel('Steps (in every 15 min)')
        % drawnow
        % 
        % subplot(2,6,2), cla
        % plot(sol_CT{2}(:),'bo-')
        % hold on 
        % plot(x_DA_int(2,:),'go-')
        % grid on
        % ylabel('H (in kWh)')
        % xlabel('Steps (in every 15 min)')
        % drawnow
        % 
        % subplot(2,6,3), cla
        % plot(sol_CT{3},'bo-')
        % hold on
        % plot(Pc_DA_int,'go-')
        % grid on
        % ylabel('Pc (in kW)')
        % xlabel('Steps (in every 15 min)')
        % drawnow
        % % 
        % subplot(2,6,4), cla
        % plot(sol_CT{4},'bo-')
        % hold on
        % plot(Pd_DA_int,'go-')
        % grid on
        % ylabel('Pd (in kW)')
        % xlabel('Steps (in every 15 min)')
        % drawnow
        % % 
        % subplot(2,6,5), cla
        % plot(sol_CT{5},'bo-')
        % hold on
        % plot(Hc_DA_int,'go-')
        % grid on
        % ylabel('Hc (in kW)')
        % xlabel('Steps (in every 15 min)')
        % drawnow
        % % 
        % subplot(2,6,6), cla
        % plot(sol_CT{6},'bo-')
        % hold on
        % plot(Hd_DA_int,'go-')
        % grid on
        % ylabel('Hd (in kW)')
        % xlabel('Steps (in every 15 min)')
        % drawnow
        % 
        % % 
        % subplot(2,6,7), cla
        % hold on
        % plot(sol_CT{13},'bo-')
        % hold on 
        % plot(Pg_DA_rshp,'go-')
        % grid on
        % ylabel('Pgrid_{CT}/Pgrid_{DA} (in kW)' )
        % xlabel('Steps (in every 15 min)')
        % drawnow
        % % 
        % subplot(2,6,8), cla
        % plot(sol_CT{8},'bo-')
        % hold on
        % plot(P_CHP_DA_int,'go-')
        % grid on
        % ylabel('P_{CHP} (in kW)')
        % xlabel('Steps (in every 15 min)')
        % drawnow
        % 
        % subplot(2,6,9), cla
        % plot(sol_CT{11},'bo-')
        % hold on
        % plot(H_EB_DA_int,'go-')
        % grid on
        % ylabel('H_{EB} (in kW)')
        % xlabel('Steps (in every 15 min)')
        % drawnow
        % 
        % subplot(2,6,10), cla
        % plot(sol_CT{12},'bo-')
        % hold on
        % plot(H_GB_DA_int,'go-')
        % grid on
        % ylabel('H_{GB} (in kW)')
        % xlabel('Steps (in every 15 min)')
        % drawnow

         % subplot(2,6,11), cla
         % plot(deltaP_CT_vec,'bo-')
         % grid on
         % ylabel('P_{WT}+P_{PV}-P_L (in kW)')
         % xlabel('Steps (in every 15 min)')
         % drawnow
       

    
 end
 toc


  

end
