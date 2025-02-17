%% This script processes the measurement data for the day ahead
%% 

clear; clc;
close all;
load ('electrical_demand_Hannover_2019');

N12 = 12;
N48 = 48;
N24 = 24;
factor = 0;%0.1;  % 10% noise level
start = (365-26)*24-12;
index = 26*N24+1+12; % start at 12:00 the day ahead

% electrical load profile
PL_DA_real = period_0(:,2)*1e-3; % from W to kW
PL_DA_real = PL_DA_real'/max(PL_DA_real)*1300; % normalize and rescale
PL_DA_real = PL_DA_real(index:end);   % we are on 27.01 at 12:00


% add noise to generate the forecast
noise1 = (rand(1,start)-0.5*ones(1,start))*2; % in [-1,1]

PL_DA = PL_DA_real + factor*PL_DA_real.*noise1;   % 10 procent
PL_DA_N24 = gen_grow_noiseDA(N24, noise1, PL_DA_real);  % account for 12h before delivery
PL_DA_N48 = gen_grow_noiseDA(N48, noise1, PL_DA_real);  


% heat load profile
load ('heat_demand_Hannover_2019');
HL_DA_real = period_0(:,2)*1e-3; % from W to kW

HL_DA_real = HL_DA_real'/max(HL_DA_real)*1500; % normalize and rescale
HL_DA_real = HL_DA_real(index:end);  % we are on 27.01 at 12:00


noise2 = (rand(1,start)-0.5*ones(1,start))*2; % in [-1,1]
HL_DA = HL_DA_real + factor*HL_DA_real.*noise2;  % 10 procent
HL_DA_N24 = gen_grow_noiseDA(N24, noise2, HL_DA_real);
HL_DA_N48 = gen_grow_noiseDA(N48, noise2, HL_DA_real);



% electricity prices
v_data_raw = readtable ('Gro_handelspreise_201901010000_201912312359_stunde');
v_DA_real = transpose(str2double(v_data_raw{:,4}))*1e-2*1e-3;   % from Mwh to Kwh
v_DA_real = v_DA_real(index:end); % we are on 27.01 at 12:00

u_DA_real = 0.9*v_DA_real;


noise3 = (rand(1,start)-0.5*ones(1,start))*2; % in [-1,1] 
% real electricity prices
v_DA = v_DA_real + factor*v_DA_real.*noise3;
v_DA_N24 = gen_grow_noiseDA(N24, noise3, v_DA_real);
v_DA_N48 = gen_grow_noiseDA(N48, noise3, v_DA_real);

u_DA = 0.9*v_DA;
u_DA_N24 = 0.9*v_DA_N24;
u_DA_N48 = 0.9*v_DA_N48;



% PV generation
load ('photovoltaic_generation_Hannover_2019');
PV_DA_real = period_0(:,2)*1; % from W to kW

PV_DA_real = PV_DA_real'*400;  % rescale
PV_DA_real = PV_DA_real(index:end); % we are on 27.01 at 12:00

noise4 = (rand(1,start)-0.5*ones(1,start))*2; % in [-1,1]
PV_DA = PV_DA_real + factor*PV_DA_real.*noise4;  % 10 procent
PV_DA_N24 = gen_grow_noiseDA(N24, noise4, PV_DA_real);
PV_DA_N48 = gen_grow_noiseDA(N48, noise4, PV_DA_real);


% Wind turbine generation
WT_data_raw =  readtable ('Wind_Turbine_Hannover_2019');
WT_DA_real = WT_data_raw{:,3};   % in Kwh
%plot(WT_DA_real)
WT_DA_real = WT_DA_real'*600;  % rescale to 1500
WT_DA_real = WT_DA_real(index:end); % we are on 27.01 at 12:00

noise5 = (rand(1,start)-0.5*ones(1,start))*2; % in [-1,1]
WT_DA = WT_DA_real; %+ factor*WT_DA_real.*noise5;   % 10 procent
WT_DA_N24 = gen_grow_noiseDA(N24, noise5, WT_DA_real);
WT_DA_N48 = gen_grow_noiseDA(N48, noise5, WT_DA_real);




%% saving


% saving noise


save("noise", "noise1", "noise2", "noise3", "noise4", "noise5")

% save DA forecast

save("DA_uncertainties.mat",...
    "PL_DA","PL_DA_N24","PL_DA_N48","HL_DA", "HL_DA_N24","HL_DA_N48", ...
    "v_DA", "v_DA_N24","v_DA_N48", "PV_DA","PV_DA_N24","PV_DA_N48",...
    "WT_DA", "WT_DA_N24", "WT_DA_N48",  "u_DA", "u_DA_N24","u_DA_N48"); 
  
save("DA_real.mat","PL_DA_real", "HL_DA_real","v_DA_real", "u_DA_real","PV_DA_real", "WT_DA_real")

%% plots

% plot the electrical power profile
figure(1)
subplot(3,1,1)
plot(PL_DA_real(1:48));
hold on 
plot(PL_DA(1:48))
hold on 
%plot(PL_DA_N24(1:48))
hold on 
%plot(PL_DA_N48(1:48))
grid on
title('DA Electrical load profile starting from 28.01.19')
xlabel('Time in hours')
ylabel('Power in kW')
xlim([1, 48])

% plot heat load 2019
subplot(3,1,2)
plot(HL_DA_real(1:48))
hold on 
plot(HL_DA(1:48))
hold on 
%plot(HL_DA_N24(1:48))
%hold on 
%plot(HL_DA_N48(1:48))
grid on
title('DA Heat load profile starting from 28.01.19')
xlabel('Time in hours')
ylabel('Power in kW')
xlim([1, 48])

% plot prices in 2019
subplot(3,1,3)
plot(v_DA_real(1:48));
hold on 
plot(v_DA(1:48))
hold on 
%plot(v_DA_N24(1:48))
% hold on 
% plot(v_DA_N48(1:48))
grid on
title('DA Prices starting from 28.01.19')
xlabel('Time in hours')
ylabel('Price in Euro/Kwh')
%xlim([1, 48])

% PV generation 
figure(2)
% plot the PV
subplot(3,1,1)
plot(PV_DA_real(1:48));
hold on 
plot(PV_DA(1:48))
hold on 
%plot(PV_DA_N24(1:48))
% hold on 
% plot(PV_DA_N48(1:48))
grid on
title('DA Photovoltaik starting from 28.01.19')
xlabel('Time in hours')
ylabel('Power in kW')
xlim([1, 48])
% 
% WT generation
subplot(3,1,2)
plot(WT_DA_real(1:48));
hold on 
plot(WT_DA(1:48))
hold on 
%plot(WT_DA_N24(1:48))
% hold on 
% plot(WT_DA_N48(1:48))
grid on
title('DA wind generation starting from 28.01.19')
xlabel('Time in hours')
ylabel('Power in kW')
xlim([1, 48])

% function  X = gen_grow_noise_DA(N, noise, real)
%     relGrowth = linspace(0.01,0.2*N/24,N); % from 1% to 20% noise at 24h
%     if char(real) == "v_DA_real"
%         relGrowth = linspace(0,0.05*N/24,N); % from 0% to 5% noise at 24h
%     end
%     l = length(real);
%     X = zeros(1,l);
%     for i=1:l/N
%         X((i-1)*N+1:i*N) = real((i-1)*N+1:i*N) +...
%         real((i-1)*N+1:i*N).*relGrowth.*noise((i-1)*N+1:i*N);
%     end
% end

function  X = gen_grow_noiseDA(N, noise, real)
    % from 0% to 20% noise at 24h
    relGrowth = linspace(0,0.2*N/24,N); 
    if char(real) == "v_DA_real"
        % from 0% to 5% noise at 24h
        relGrowth = linspace(0,0.05*N/24,N); 
    end
    i = 0;
    for j= 1:24:length(real)-N+1
        i = i+1;
        X(i,:) = real(j:j+N-1) + real(j:j+N-1).*relGrowth.*noise(j:j+N-1);
    end
end
