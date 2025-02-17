 close all   

%% load measurement data
% 
% load ('...\data\8760h-01p\DA_uncertainties.mat');
% load ('...\data\8760h-01p\CT_uncertainties.mat');
% load ('...\data\8760h-01p\CT_real.mat');
% load ('...\data\8760h-01p\DA_real.mat');

   
  PL_forecast = []; 
  HL_forecast = [];
  WT_forecast = [];
  PV_forecast = [];
  v_forecast = [];


  PL_forecast = PL_DA_N48(1,1:48);
  HL_forecast = HL_DA_N48(1,1:48);
  WT_forecast = WT_DA_N48(1,1:48);
  PV_forecast = PV_DA_N48(1,1:48);
  v_forecast = v_DA_N48(1,1:48);
  u_forecast = u_DA_N48(1,1:48);

  
% plot the electrical power profile
figure(6)
subplot(2,2,1)
plot(PL_DA_real(1:48));
hold on 
plot(PL_forecast)
hold on 
grid on
xlabel('Time in hours')
ylabel('PL in kW')
xlim([1, 48])

% plot heat load 2019
subplot(2,2,2)
plot(HL_DA_real(1:48))
hold on 
plot(HL_forecast)
hold on 
grid on
xlabel('Time in hours')
ylabel('HL in kW')
xlim([1, 48])
% 
% plot prices in 2019


% plot the PV
subplot(2,2,3)
plot(PV_DA_real);
hold on 
plot(PV_forecast)
hold on 
grid on
xlabel('Time in hours')
ylabel('PV in kW')
xlim([1, 48])
% 
% WT generation
subplot(2,2,4)
plot(WT_DA_real);
hold on 
plot(WT_forecast)
hold on 
grid on
%title('DA wind generation starting from 28.01.19')
xlabel('Time in hours')
ylabel('WT in kW')
xlim([1, 48])

 % cleanfigure;
 % matlab2tikz('gen_load_DA.tex', 'parseStrings',false)

figure(7)
subplot(3,1,1)
plot(v_DA_real(1:48));
hold on 
plot(v_forecast)
hold on 
grid on
xlabel('Time in hours')
ylabel('Price in \$/Kwh')
xlim([1, 48])

subplot(3,1,2)
plot(u_DA_real(1:48));
hold on 
plot(u_forecast)
hold on 
grid on
xlabel('Time in hours')
ylabel('Price in \$/Kwh')
xlim([1, 48])

cleanfigure;
matlab2tikz('prices_DA.tex', 'parseStrings',false)

