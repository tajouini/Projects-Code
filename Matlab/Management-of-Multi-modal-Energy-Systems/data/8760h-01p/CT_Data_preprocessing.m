%% This script processes the measurement data for the continuous trading
%% 

function CT_Data_preprocessing (factor)

% if exist ("CT_real.mat", 'file') && exist ("CT_uncertainties.mat", 'file')
%   return;
% end;

%clear; clc;
close all;

N1 = 1;
N24 = 24;
N4 = 4;
N96 = 96;
N48 = 12*4;
%factor = 0.1;

%% load raw data for continuous trading
% electrical load demand
load ('DA_uncertainties.mat');
load ('DA_real.mat');
% load noise sequences
load ('noise');
% electricity prices
if exist ('OCTAVE_HOME','builtin') == 5
    v_data_CT = csv2cell ('Gro_handelspreise_201901010000_201912312359_Viertelstunde.csv',1,';');
    v_CT_real = transpose(cell2mat(v_data_CT(:,4)));   
else
    v_data_CT = readtable ('Gro_handelspreise_201901010000_201912312359_Viertelstunde');
    v_CT_real = transpose(v_data_CT{:,4});   
end
% start as in DA: on 27.01 at 11:45
index = (26*N24+1+11)*4+1;
% interpolated v_DA_real
v_CT_real = v_CT_real(index:end)*1e-3; % from Mwh to Kwh;  
%% linear interpolation of the raw data 
% CT electrical load profile
% interpolating between the hourly power values to get to every 15 min
PL_CT_real=zeros(1,4*length(PL_DA_real));
HL_CT_real=zeros(1,4*length(HL_DA_real));
PV_CT_real = zeros(1,4*length(PV_DA_real));
WT_CT_real=zeros(1,4*length(WT_DA_real));

noise1_CT = zeros(1,4*length(noise1));
noise2_CT = zeros(1,4*length(noise2));
noise3_CT = zeros(1,4*length(noise3));
noise4_CT = zeros(1,4*length(noise4));
noise5_CT = zeros(1,4*length(noise5));

for i=1:1:length(PL_DA_real)-1
  d = 1*(PL_DA_real(i+1)-PL_DA_real(i))/4;
  PL_CT_real(4*(i-1)+1) = PL_DA_real(i);
  PL_CT_real(4*(i-1)+2) = PL_DA_real(i) + d;
  PL_CT_real(4*(i-1)+3) = PL_DA_real(i) + d*2;
  PL_CT_real(4*i) = PL_DA_real(i) +d*3;
  
  d = 1*(noise1(i+1)-noise1(i))/4;
  noise1_CT(4*(i-1)+1) = noise1(i);
  noise1_CT(4*(i-1)+2) = noise1(i) + d;
  noise1_CT(4*(i-1)+3) = noise1(i) + d*2;
  noise1_CT(4*i) = noise1(i) +d*3;

  d = 1*(HL_DA_real(i+1)-HL_DA_real(i))/4;
  HL_CT_real(4*(i-1)+1) = HL_DA_real(i);
  HL_CT_real(4*(i-1)+2) = HL_DA_real(i) + d;
  HL_CT_real(4*(i-1)+3) = HL_DA_real(i) + d*2;
  HL_CT_real(4*i) = HL_DA_real(i) +d*3;
  
  d = 1*(noise2(i+1)-noise2(i))/4;
  noise2_CT(4*(i-1)+1) = noise2(i);
  noise2_CT(4*(i-1)+2) = noise2(i) + d;
  noise2_CT(4*(i-1)+3) = noise2(i) + d*2;
  noise2_CT(4*i) = noise2(i)+ d*3;  
  
  d = 1*(noise3(i+1)-noise3(i))/4;
  noise3_CT(4*(i-1)+1) = noise3(i);
  noise3_CT(4*(i-1)+2) = noise3(i) + d;
  noise3_CT(4*(i-1)+3) = noise3(i) + d*2;
  noise3_CT(4*i) = noise3(i)+ d*3;
  
  d = 1*(PV_DA_real(i+1)-PV_DA_real(i))/4;
  PV_CT_real(4*(i-1)+1) = PV_DA_real(i);
  PV_CT_real(4*(i-1)+2) = PV_DA_real(i) + d;
  PV_CT_real(4*(i-1)+3) = PV_DA_real(i) + d*2;
  PV_CT_real(4*i) = PV_DA_real(i) +d*3;
  
  d = 1*(noise4(i+1)-noise4(i))/4;
  noise4_CT(4*(i-1)+1) = noise4(i);
  noise4_CT(4*(i-1)+2) = noise4(i) + d;
  noise4_CT(4*(i-1)+3) = noise4(i) + d*2;
  noise4_CT(4*i) = noise4(i)+ d*3;
  
  d = 1*(WT_DA_real(i+1)-WT_DA_real(i))/4;
  WT_CT_real(4*(i-1)+1) = WT_DA_real(i);
  WT_CT_real(4*(i-1)+2) = WT_DA_real(i) + d;
  WT_CT_real(4*(i-1)+3) = WT_DA_real(i) + d*2;
  WT_CT_real(4*i) = WT_DA_real(i) +d*3;
  
  d = 1*(noise5(i+1)-noise5(i))/4;
  noise5_CT(4*(i-1)+1) = noise5(i);
  noise5_CT(4*(i-1)+2) = noise5(i) + d;
  noise5_CT(4*(i-1)+3) = noise5(i) + d*2;
  noise5_CT(4*i) = noise5(i)+ d*3;
end

for i = 0:3 
  PL_CT_real(end-i) = PL_DA_real(end);
  HL_CT_real(end-i) = HL_DA_real(end);
  PV_CT_real(end-i) = PV_DA_real(end);
  WT_CT_real(end-i) = WT_DA_real(end);
end

% start on 27.01.23 at 23:45
PL_CT_real = PL_CT_real(12*4:end);
HL_CT_real = HL_CT_real(12*4:end);
PV_CT_real = PV_CT_real(12*4:end);
WT_CT_real = WT_CT_real(12*4:end);
v_CT_real = v_CT_real(12*4:end);  
v_CT_real = 1.2*v_CT_real; %1.2*v_CT_real;  % 1.2*v_DA

u_CT_real = 0.8*0.9/1.2*v_CT_real; %0.8*0.9/1.2*v_CT_real;   

noise1_CT = noise1_CT(12*4:end);
noise2_CT = noise2_CT(12*4:end);
noise3_CT = noise3_CT(12*4:end);
noise4_CT = noise4_CT(12*4:end);
noise5_CT = noise5_CT(12*4:end);

%% CT forcast data 
PL_CT = PL_CT_real; 

PL_CT_N4 = gen_grow_noise_CT(N4+N1, noise1_CT, PL_CT_real, factor);
PL_CT_N96 = gen_grow_noise_CT(N96+N1, noise1_CT, PL_CT_real, factor);
PL_CT_N48 = gen_grow_noise_CT(N48+N1, noise1_CT, PL_CT_real, factor);


HL_CT = HL_CT_real; 
HL_CT_N4 = gen_grow_noise_CT(N4+N1, noise2_CT, HL_CT_real, factor);
HL_CT_N96 = gen_grow_noise_CT(N96+N1, noise2_CT, HL_CT_real, factor);
HL_CT_N48 = gen_grow_noise_CT(N48+N1, noise2_CT, HL_CT_real, factor);

v_CT = v_CT_real; 
v_CT_N4 = gen_grow_noise_CT(N4+N1, noise3_CT, v_CT_real, factor);
v_CT_N48 = gen_grow_noise_CT(N48+N1, noise3_CT, v_CT_real, factor);
v_CT_N96 = gen_grow_noise_CT(N96+N1, noise3_CT, v_CT_real, factor);

u_CT = 0.8*0.9/1.2*v_CT;
u_CT_N4 = 0.8*0.9/1.2*v_CT_N4;
u_CT_N48 = 0.8*0.9/1.2*v_CT_N48;
u_CT_N96 = 0.8*0.9/1.2*v_CT_N96;


PV_CT = PV_CT_real; 
PV_CT_N4 = gen_grow_noise_CT(N4+N1, noise4_CT, PV_CT_real, factor);
PV_CT_N48 = gen_grow_noise_CT(N48+N1, noise4_CT, PV_CT_real, factor);
PV_CT_N96 = gen_grow_noise_CT(N96+N1, noise4_CT, PV_CT_real, factor);

WT_CT = WT_CT_real; 
WT_CT_N4 = gen_grow_noise_CT(N4+N1, noise5_CT, WT_CT_real, factor);
WT_CT_N48 = gen_grow_noise_CT(N48+N1, noise5_CT, WT_CT_real, factor);
WT_CT_N96 = gen_grow_noise_CT(N96+N1, noise5_CT, WT_CT_real, factor);

%% save data

save("CT_uncertainties.mat",...
    "PL_CT","PL_CT_N4","PL_CT_N48", "PL_CT_N96","HL_CT", "HL_CT_N4","HL_CT_N48", "HL_CT_N96", ...
    "v_CT", "v_CT_N4","v_CT_N48","v_CT_N96", "PV_CT","PV_CT_N4","PV_CT_N48", "PV_CT_N96",...
    "WT_CT", "WT_CT_N4", "WT_CT_N48", "WT_CT_N96",  "u_CT", "u_CT_N4","u_CT_N48", "u_CT_N96", '-v6'); 

save("CT_real.mat","PL_CT_real", "HL_CT_real","v_CT_real","u_CT_real", "PV_CT_real", "WT_CT_real", '-v6')

%% Plots
figure(1)
% time in min
t= 1:length(PL_CT(1:48*4)); 
subplot(2,1,1)
plot(PL_CT_real(1:48*4));
hold on 
plot(t,PL_CT(1:48*4));
hold on 
%plot(t,PL_CT_N4(1:48*4));
hold on 
%plot(t,PL_CT_N96(1:48*4));
grid on
legend('Pload', 'Hload')
title('CT Load profiles')
xlabel('Step (every 15 min)')
ylabel('Power in kW')

subplot(2,1,2)
plot(t,HL_CT_real(1:48*4));
hold on 
plot(t,HL_CT(1:48*4));
grid on
legend('Hload real', 'Hload')
title('CT Load profiles')
xlabel('Step (every 15 min)')
ylabel('Power in kW')
% xlim([1, length(PL_CT)])


figure(2)
plot(t,v_CT_real(1:48*4));
hold on 
plot(t,v_CT(1:48*4));
grid on
title('CT Purchasing Price profiles')
xlabel('Step (every 15 min)')
ylabel('Power in kW')

figure(3)
subplot(2,1,1)
plot(PV_CT_real(1:48*4));
hold on 
plot(PV_CT(1:48*4));
grid on
legend('PV real', 'PV')
title('PV generation of CT profiles')
xlabel('Step (every 15 min)')
ylabel('Price in Euro')

subplot(2,1,2)
plot(WT_CT_real(1:48*4));
hold on 
plot(WT_CT(1:48*4));
grid on
legend('WT real', 'WT')
title('WT generation of CT profiles')
xlabel('Step (every 15 min)')
ylabel('Price in Euro')

end


function  X = gen_grow_noise_CT(N, noise, real, factor)
    relGrowth = linspace(0,factor*N/(24*4),N); % from 0% to 20% noise at-- 24h

    if strcmp (inputname(3), "v_DA_real")
        relGrowth = linspace(0,0.05*N/(24*4),N); % from 0% to 5% noise at 24h
    end
    for j= 1:length(real)-N+1
        X(j,:) = real(j:j+N-1) + real(j:j+N-1).*relGrowth.*noise(j:j+N-1);
    end
end
