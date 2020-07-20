clc
clear all
%% Data Inputs 
%ip = importdata('F:\Earthquake data\Magnitude 7\All Mag 7 vt\TW_vt.txt'); % Data file path
Acc_EW = importdata('D:\Fun Progs\MATLAB 3D  and STA-LTA data plus code\data plus code\Nepal\ADIB.HHE.dat');
Acc_NS = importdata('D:\Fun Progs\MATLAB 3D  and STA-LTA data plus code\data plus code\Nepal\ADIB.HHN.dat');
Acc_ver = importdata('D:\Fun Progs\MATLAB 3D  and STA-LTA data plus code\data plus code\Nepal\ADIB.HHZ.dat');
Fs = 200; %sampling frequency
%% Signal Pre-Processing
%Filter Design
digfilt = designfilt('lowpassiir', 'PassbandFrequency', 20, 'StopbandFrequency', 25, 'PassbandRipple', 1, 'StopbandAttenuation', 60, 'SampleRate', 200);
% Filtering Data
Acc_EW_filt = filter(digfilt,Acc_EW);
Acc_NS_filt = filter(digfilt,Acc_NS);
Acc_ver_filt = filter(digfilt,Acc_ver);
Fhp = 0.8; % high pass filter cutofff frequency 
[b1,a1] = butter(3,Fhp/Fs,'high'); %3rd order high pass butterworth filter
fildat = filter(b1,a1,Acc_ver); % filtered acceleration data
vel = cumtrapz(fildat)./Fs; % Integrating acceleration data for velocity
[b2,a2] = butter(3,Fhp/Fs,'high'); %3rd order high pass butterworth filter
fildat1 = filter(b2,a2,vel); % filtered velocity data
dis = cumtrapz(fildat1)./Fs;  % Integrating velocity data for displacement
peakToPeakRange = max(fildat) - min(fildat);
dt = 1/Fs; %sampling time
nt = length(fildat); % length of the input signal
time = (1:nt).*dt; % time duration of the input signal
%% Removing the delay introduced by lowpassiir filter
Acc_EW_dlycompensated = filtfilt(digfilt,Acc_EW);
Acc_NS_dlycompensated = filtfilt(digfilt,Acc_NS);
Acc_ver_dlycompensated = filtfilt(digfilt,Acc_ver);
 %% Finding the SNR
%  noise_ratio = snr(fildat) %returns the SNR in dBc of a real sinusoidal input signal, x, sampled at a rate fs. 
%  % The computation excludes the power contained in the lowest n harmonics, including the fundamental. 
%  % The default value of fs is 1. The default value of n is 6.
%  %fildat( .*) = 10^(3/20); %To increase the power of x by 3 dB:

%% plotting frequency spectrum
% fftsig = fft(fildat); %Taking fourier transform
% fftSig = fftshift(fftsig); %apply fftshift to put it in the form we are used to
% f = [-nt/2:nt/2-1]/nt; %Next, calculate the frequency axis, which is defined by the sampling rate
% figure;
% plot(f, abs(fftSig));
% title('magnitude FFT of sine');
% xlabel('Frequency (Hz)');
% ylabel('magnitude');
%% STA-LTA Algorithm gor P-Wave detection
stw = 1;    %short time window length
ltw = 60;   %long time window length
thresh = 4 ; % Threshold
thresh1 = 3;
%t = 1;      
nl = fix(ltw / dt); %no. of data points in the long time window
ns = fix(stw / dt); %no. of data points in the short time window 
nt = length(fildat); 
sra = zeros(1, nt);
for k = nl+1:nt
    sta(k,1) = (1/ns)* trapz(abs(fildat(k-ns:k)));
    lta(k,1) = (1/nl)* trapz(abs(fildat(k-nl:k)));
 end
for l = nl+1: nt 
    sra(l) = sta(l)/lta(l);
end
  itm = find(sra > thresh);
    if ~isempty(itm)
      itmax = itm(1);
    end
    tp = itmax*dt; % P-wave arriving time 
    fprintf('P-Wave detection time for threshold 4 =  %f second\n', tp);
    
  itm1 = find(sra > thresh1);
    if ~isempty(itm1)
      itmax1 = itm1(1);
    end
    tp1 = itmax1*dt; % P-wave arriving time 
    fprintf('P-Wave detection time for threshold 3 = %f second\n', tp1);
%% S-wave arrival time
pkHts = 0.72; % 10 percent
[pk2,t22] = findpeaks(Acc_NS_dlycompensated,Fs,'MinPeakHeight',pkHts*max(Acc_ver_dlycompensated),'Npeaks',1);
[pk3,t33] = findpeaks(Acc_EW_dlycompensated,Fs,'MinPeakHeight',pkHts*max(Acc_ver_dlycompensated),'Npeaks',1);

display(sprintf('S-wave found on EW component at %f seconds and on NS componet at %f seconds,', t33,t22));

if(t22<t33)
    display('S-wave detected first on North-South component');
else
    display('S-wave detected first on East-West component');
end
ts = min(t22,t33);
line([ts,ts],[min(get(gca,'Ylim'))],'linestyle','--','linewidth',2,'color','red');
    
%% Tauc , Pd and Magnitude calculations
vel_sq = vel.^2;
dis_sq = dis.^2;
r1 = trapz(vel_sq((itmax):(itmax+600)));
r2 = trapz(dis_sq((itmax):(itmax+600)));
r = r1/r2;
tauc = 2*pi/sqrt(r);
pd = max(dis((itmax):(itmax+600)));
% mag_tauc = (log(tauc) + 3.45)/0.47 %Coefficients varies from region to region
% mag_pd = (0.873*((log(pd)+6.3)/0.513))+4.74 %Coefficients varies from region to region

%% Distance of earthquake from the seismometer
dist = (ts-tp)*8;
display(sprintf('Earthquake is estimated to be %f kilometers from the seismometer',dist))

%% Acceleration Plot
figure(1);
subplot(3,1,1)
plot(time,fildat,[tp tp],ylim,'r','LineWidth',1)
%plot(time,fildat)
title('Acceleration Data');
xlabel('Time (Sec)');
ylabel('Acceleration (cm/sec^2)');
grid on 
grid minor
%% Velocity Plot
subplot(3,1,2)
plot(time,vel)
title('Velocity Data');
xlabel('Time (Sec)');
ylabel('Velocity (cm/sec)');
grid on 
grid minor
%% Displacement Plot
subplot(3,1,3)
plot(time,dis)
title('Displacement Data');
xlabel('Time (Sec)');
ylabel('Displacement (cm)');
grid on 
grid minor

%% Plotting Spectogram of Original Signal and detecting the P-wave first arrival
figure(2)
box on
hold on
subplot(3,1,1)
plot(time,fildat,[tp tp],ylim,'r','LineWidth',2)
hold on
plot(time,fildat,[tp1 tp1],ylim,'g','LineWidth',2)
%line([tp1 tp1],[0,100],'Color','green','LineWidth',2);
%title('Acceleration Data');
xlabel('Time (Sec)');
ylabel('Acceleration (cm/sec^2)');
grid on 
grid minor
axis tight

box on
s = spectrogram(abs(fildat),256,250,256,200,'yaxis');
subplot(3,1,2)
%title('Spectrogram of Acceleration Data');
spectrogram(abs(fildat),256,250,256,200,'yaxis')
tp_in_min = tp/60;
tp_in_min1 = tp1/60;
%line([tp_in_min tp_in_min],[0,100],'Color','red','LineWidth',2);
line([tp_in_min1 tp_in_min1],[0,100],'Color','green','LineWidth',2);
grid on 
grid minor
axis tight;

box on
subplot(3,1,3)
thresh_spec = spectrogram(abs(fildat),256,250,256,200,'MinThreshold',-50,'yaxis');
thresh_spec1 = abs(thresh_spec);
%title('Spectrogram of Acceleration Data with -50dB Threshold');
spectrogram(abs(fildat),256,250,256,200,'MinThreshold',-50,'yaxis')
tp_in_min = tp/60;
tp_in_min1 = tp1/60;
%line([tp_in_min tp_in_min],[0,100],'Color','red','LineWidth',2);
line([tp_in_min1 tp_in_min1],[0,100],'Color','green','LineWidth',2);

grid on 
grid minor
axis tight;

%print('-dpdf','-r600','Temp.pdf')

%% Vizualizing the displacement of the seismometer

% Using an integrator to compute velocity from acceleration data.
% Integrator is filter with TF = 1/(1-Z^-1)
Nr = 1;
Dr = [1,-1];
% Detrending acceleration data it to remove drift,Integrating and dividing with sampling frequency to get velocity data
Ver_acc_dt = detrend(Acc_ver_dlycompensated);
NS_acc_dt = detrend(Acc_NS_dlycompensated);
EW_acc_dt = detrend(Acc_EW_dlycompensated);

vel_ver_nodrift = filter(Nr,Dr,Ver_acc_dt)/Fs; 
vel_NS_nodrift = filter(Nr,Dr,NS_acc_dt)/Fs;
vel_EW_nodrift = filter(Nr,Dr,EW_acc_dt)/Fs;

% Integrating velocity data and dividing with sampling frequency to get displacement data
dis_ver = filter(Nr,Dr,vel_ver_nodrift)/Fs;
dis_NS = filter(Nr,Dr,vel_NS_nodrift)/Fs;
dis_EW = filter(Nr,Dr,vel_EW_nodrift)/Fs;

% Plot the positions using 3D plot to visualize seismometer movements
figure(3);
plot3(dis_NS,dis_EW,dis_ver);
grid on; view([-45,30]);
xlabel('N-S Direction in cm');
ylabel('E-W Direction in cm');
zlabel('Vertical Direction in cm');
title('Displacement of the seismometer in 3D');
set(gcf,'Name','Seismometer Trajectory');
set(gcf,'Units','Normalized');