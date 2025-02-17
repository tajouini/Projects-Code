%% run the hiearchichal MPC scheme with different horizon lengths



% generate a new sequence
% cd('../data/8760h-01p')
% DA_Data_preprocessing
% CT_Data_preprocessing
% cd('../../code')
clear
close all
addpath('c:\gurobi1002\win64\matlab')

% no noise
days = 4;


% noisy case  
 noise = 1;
% 
 [index_DA, index_mon_DA, index_mon_CT, index_CT, index_real_t, C_input,flag_DA, flag_CT] = ...
  MPC_hiearchical_approach_third_layer(24,24,1/4,days,noise);
% 
%  
 % [index_DA, index_mon_DA, index_mon_CT, index_CT, index_real_t, flag_DA, flag_CT] = ...
 % MPC_hiearchical_approach_third_layer(48,24,1/4,days,noise);
% 

% [index_DA_48_24_noise, index_CT_48_24_noise, index_mon_DA_48_24_noise, index_mon_CT_48_24_noise] = ...
% MPC_hiearchical_approach_fix_u1(48,24,days,noise);

%[index_DA_24_12, index_CT_24_12, index_mon_DA_24_12, index_mon_CT_24_12] = ...
%MPC_hiearchical_approach_fix_u1(24,12,days,noise);

