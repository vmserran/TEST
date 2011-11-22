%Program to read pressure, resistance and temperature data from BioTAC's
%Change of directory
cd('C:\Users\VICKY\Documents\Haptic Library\BIOTAC SIGNALS USING CHEETAH AND C++\BIOTACS AT REST\TRIAL 1');
clear all; close all; clc;

%% AC Pressure

% AC pressure raw data
Pac_bt1_raw = textread('ac_pressure_biotac_1.txt');
Pac_bt2_raw = textread('ac_pressure_biotac_2.txt');
Pac_bt3_raw = textread('ac_pressure_biotac_3.txt');

% retrieving the AC Pressure data-Biotac 1
Pac_sampling_freq = 2200;
num_of_Pac_readings1 = length(Pac_bt1_raw);
Pac_time_end1 = num_of_Pac_readings1/Pac_sampling_freq;
Pac_time1 = linspace(0, Pac_time_end1-(1/Pac_sampling_freq), num_of_Pac_readings1);%Reading AC Pressure from all 3 Biotacs-raw data

% retrieving the AC Pressure data-Biotac 2
num_of_Pac_readings2 = length(Pac_bt2_raw);
Pac_time_end2 = num_of_Pac_readings2/Pac_sampling_freq;
Pac_time2 = linspace(0, Pac_time_end2-(1/Pac_sampling_freq), num_of_Pac_readings2);

% retrieving the AC Pressure data-Biotac 3
num_of_Pac_readings3 = length(Pac_bt3_raw);
Pac_time_end3 = num_of_Pac_readings3/Pac_sampling_freq;
Pac_time3 = linspace(0, Pac_time_end3-(1/Pac_sampling_freq), num_of_Pac_readings3);

%Converting raw data into pressure in upsi units
Pac_bt1=(Pac_bt1_raw-mean(Pac_bt1_raw(1:100))).*54.211762;
Pac_bt2=(Pac_bt2_raw-mean(Pac_bt2_raw(1:100))).*54.211762;
Pac_bt3=(Pac_bt3_raw-mean(Pac_bt3_raw(1:100))).*54.211762;

% Plots
figure (1); hold on; grid on;
subplot(3,1,1);plot (Pac_time1,Pac_bt1,'b');title('AC Pressure for BT1 using C++');xlabel('time (s)');ylabel('P(upsi)');
subplot(3,1,2);plot (Pac_time2,Pac_bt2,'r');title('AC Pressure for BT2 using C++');xlabel('time (s)');ylabel('P(upsi)');
subplot(3,1,3);plot (Pac_time3,Pac_bt3,'g');title('AC Pressure for BT3 using C++');xlabel('time (s)');ylabel('P(upsi)');
%legend ('BT1','BT2','BT3')


%% Impedances

%Impedance values, TDC, TAC, PDC from Biotac 1
Imp_raw_BT1 = textread('dc_imp_values_1.txt');
numrow_BT1 = 22;
size(Imp_raw_BT1,1);
numcol_BT1 =(size(Imp_raw_BT1,1)/numrow_BT1);

%Reshaping the information obtained in textfile
if numcol_BT1=='round'
    numcol_BT1=(size(Imp_raw_BT1,1)/numrow_BT1);
    Imp_reshape_BT1 = reshape(Imp_raw_BT1(1:numrow_BT1*numcol_BT1), numrow_BT1, numcol_BT1);
else
    numcol_BT1=floor(size(Imp_raw_BT1,1)/numrow_BT1);
    Imp_reshape_BT1 = reshape(Imp_raw_BT1(1:numrow_BT1*numcol_BT1), numrow_BT1, numcol_BT1);
end

Imp_BT1 = Imp_reshape_BT1(1:19,:);

%Converting raw data into resistance in k-ohms
R_BT1=(40950-10*Imp_BT1)./Imp_BT1;
 
%Impedance values, TDC, TAC, PDC from Biotac 2
Imp_raw_BT2 = textread('dc_imp_values_2.txt');
numrow_BT2 = 22;
size(Imp_raw_BT2,1);
numcol_BT2 =(size(Imp_raw_BT2,1)/numrow_BT2);

%Reshaping the information obtained in textfile
if numcol_BT2=='round'
    numcol_BT2=(size(Imp_raw_BT2,1)/numrow_BT2);
    Imp_reshape_BT2 = reshape(Imp_raw_BT2(1:numrow_BT2*numcol_BT2), numrow_BT2, numcol_BT2);
else
    numcol_BT2=floor(size(Imp_raw_BT2,1)/numrow_BT2);
    Imp_reshape_BT2 = reshape(Imp_raw_BT2(1:numrow_BT2*numcol_BT2), numrow_BT2, numcol_BT2);
end

Imp_BT2 = Imp_reshape_BT2(1:19,:);

%Converting raw data into resistance in k-ohms
R_BT2=(40950-10*Imp_BT2)./Imp_BT2;

%Impedance values, TDC, TAC, PDC from Biotac 3
Imp_raw_BT3 = textread('dc_imp_values_3.txt');
numrow_BT3 = 22;
size(Imp_raw_BT3,1);
numcol_BT3 =(size(Imp_raw_BT3,1)/numrow_BT3);

%Reshaping the information obtained in textfile
if numcol_BT3=='round'
    numcol_BT3=(size(Imp_raw_BT3,1)/numrow_BT3);
    Imp_reshape_BT3 = reshape(Imp_raw_BT3(1:numrow_BT3*numcol_BT3), numrow_BT3, numcol_BT3);
else
    numcol_BT3=floor(size(Imp_raw_BT3,1)/numrow_BT3);
    Imp_reshape_BT3 = reshape(Imp_raw_BT3(1:numrow_BT3*numcol_BT3), numrow_BT3, numcol_BT3);
end

Imp_BT3 = Imp_reshape_BT3(1:19,:);

%Converting raw data into resistance in k-ohms
R_BT3=(40950-10*Imp_BT3)./Imp_BT3;

%Time scale
electrodes_sampling_freq = 100;
electrodes_time1 = 0 : 1/electrodes_sampling_freq : ...
                        (length(Imp_BT1)-1)/electrodes_sampling_freq;
electrodes_time2 = 0 : 1/electrodes_sampling_freq : ...
                        (length(Imp_BT2)-1)/electrodes_sampling_freq;
electrodes_time3 = 0 : 1/electrodes_sampling_freq : ...
                        (length(Imp_BT3)-1)/electrodes_sampling_freq;
                    
%Plots
figure (2)
subplot (3,1,1); plot(electrodes_time1,R_BT1');title('Impedance values for Biotac 1 using C++');xlabel('data points');ylabel('Resistance(Kilo-ohm)')
subplot (3,1,2); plot(electrodes_time2,R_BT2');title('Impedance values for Biotac 2 using C++');xlabel('data points');ylabel('Resistance(Kilo-ohm)')
subplot (3,1,3); plot(electrodes_time3,R_BT3');title('Impedance values for Biotac 3 using C++');xlabel('data points');ylabel('Resistance(Kilo-ohm)')
legend ('data1','data2','data3','data4','data5','data6','data7','data8','data9','data10','data11','data12','data13','data14','data15','data16','data17','data18','data19');

%% DC Pressure     
PDC_raw_BT1 = Imp_reshape_BT1(20,:);
PDC_raw_BT2 = Imp_reshape_BT2(20,:);
PDC_raw_BT3 = Imp_reshape_BT3(20,:);

%Converting raw data into  psi units
PDC_BT1=(PDC_raw_BT1-mean(PDC_raw_BT1(1:100)))*0.005372;
PDC_BT2=(PDC_raw_BT2-mean(PDC_raw_BT2(1:100)))*0.005372;
PDC_BT3=(PDC_raw_BT3-mean(PDC_raw_BT3(1:100)))*0.005372;

%Time scale
Pdc_sampling_freq = 100;
Pdc_time1 = 0 : 1/Pdc_sampling_freq : ...
                        (length(PDC_BT1)-1)/Pdc_sampling_freq;
Pdc_time2 = 0 : 1/Pdc_sampling_freq : ...
                        (length(PDC_BT3)-1)/Pdc_sampling_freq;
Pdc_time3 = 0 : 1/Pdc_sampling_freq : ...
                        (length(PDC_BT3)-1)/Pdc_sampling_freq;
                    
%Plots
figure (3); hold on; grid on
plot (Pdc_time1,PDC_BT1,'b');title('DC Pressure using C++');xlabel('time (s)');ylabel('P(psi)');
plot (Pdc_time2,PDC_BT2, 'r');
plot (Pdc_time3,PDC_BT3, 'g');
ylim([-0.01 0.015])
legend ('BT1','BT2','BT3')

 
%% AC Temperature 
TAC_raw_BT1 = Imp_reshape_BT1(21,:);
TAC_raw_BT2 = Imp_reshape_BT2(21,:);
TAC_raw_BT3 = Imp_reshape_BT3(21,:);

%Converting raw data into Celsius units
TAC_BT1=25.59./(TAC_raw_BT1-mean(TAC_raw_BT1(1:100))).^2;
TAC_BT2=25.59./(TAC_raw_BT2-mean(TAC_raw_BT2(1:100))).^2;
TAC_BT3=25.59./(TAC_raw_BT3-mean(TAC_raw_BT3(1:100))).^2;

%Plots
figure (4)
subplot (3,1,1);plot (TAC_BT1);title('BioTac 1-AC Temperature');xlabel('time (s)');ylabel('T(oC)');
subplot (3,1,2);plot (TAC_BT2);title('BioTac 2-AC Temperature');xlabel('time (s)');ylabel('T(oC)');
subplot (3,1,3);plot (TAC_BT3);title('BioTac 3-AC Temperature');xlabel('time (s)');ylabel('T(oC)');

figure (5)
hold on
plot (TAC_BT1,'b');plot (TAC_BT2,'g');plot (TAC_BT3,'r');
title('AC Temperature');xlabel('time (s)');ylabel('T(oC)');
legend ('BT1','BT2','BT3')

%% DC Temperature 
TDC_raw_BT1 = Imp_reshape_BT1(22,:);
TDC_raw_BT2 = Imp_reshape_BT2(22,:);
TDC_raw_BT3 = Imp_reshape_BT3(22,:);

%Converting raw data into Celsius units
TDC_BT1=-25.59*log(TDC_raw_BT1)+209.46;
TDC_BT2=-25.59*log(TDC_raw_BT2)+209.46;
TDC_BT3=-25.59*log(TDC_raw_BT3)+209.46;

%Plots
figure (6)
subplot (3,1,1);plot (TDC_BT1);title('BioTac 1-DC Temperature');xlabel('time (s)');ylabel('T(oC)');
subplot (3,1,2);plot (TDC_BT2);title('BioTac 2-DC Temperature');xlabel('time (s)');ylabel('T(oC)');
subplot (3,1,3);plot (TDC_BT3);title('BioTac 3-DC Temperature');xlabel('time (s)');ylabel('T(oC)');

figure (7)
hold on
plot (TDC_BT1,'b');plot (TDC_BT2,'g');plot (TDC_BT3,'r');
title('DC Temperature');xlabel('time (s)');ylabel('T(oC)');
legend ('BT1','BT2','BT3')
