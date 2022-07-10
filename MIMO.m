clc
clear
close all

c = physconst('lightspeed');
BW = 2e9;
MagTX1=load('MagsTX1.tab');
PhaseTX1=load('PhasesTX1.tab');

MagTX2=load('MagsTX2.tab');
PhaseTX2=load('PhasesTX2.tab');

MagTX3=load('MagsTX3.tab');
PhaseTX3=load('PhasesTX3.tab');

MagTX4=load('MagsTX4.tab');
PhaseTX4=load('PhasesTX4.tab');

TX1RX123=MagTX1(:,2:end).*exp(j*PhaseTX1(:,2:end));
TX2RX123=MagTX2(:,2:end).*exp(j*PhaseTX2(:,2:end));
TX3RX123=MagTX3(:,2:end).*exp(j*PhaseTX3(:,2:end));
TX4RX123=MagTX4(:,2:end).*exp(j*PhaseTX4(:,2:end));

TX1RX1=TX1RX123(:,1);
TX1RX2=TX1RX123(:,2);
TX1RX3=TX1RX123(:,3);

TX2RX1=TX2RX123(:,1);
TX2RX2=TX2RX123(:,2);
TX2RX3=TX2RX123(:,3);

TX3RX1=TX3RX123(:,1);
TX3RX2=TX3RX123(:,2);
TX3RX3=TX3RX123(:,3);

TX4RX1=TX4RX123(:,1);
TX4RX2=TX4RX123(:,2);
TX4RX3=TX4RX123(:,3);

N=length(TX1RX123);

Delta=c/(2*BW);
R_max=(N-1)*Delta;

r=linspace(0,R_max,N*10);
rA=linspace(-pi/2,pi/2,12);

T1R1 = ifft(TX1RX1,N*10);
T1R2 = ifft(TX1RX2,N*10);
T1R3 = ifft(TX1RX3,N*10);

T2R1 = ifft(TX2RX1,N*10);
T2R2 = ifft(TX2RX2,N*10);
T2R3 = ifft(TX2RX3,N*10);

T3R1 = ifft(TX3RX1,N*10);
T3R2 = ifft(TX3RX2,N*10);
T3R3 = ifft(TX3RX3,N*10);

T4R1 = ifft(TX4RX1,N*10);
T4R2 = ifft(TX4RX2,N*10);
T4R3 = ifft(TX4RX3,N*10);

figure 
plot(r,10*log10(abs(T1R1)));
xlabel("Range [m]");
ylabel("P [dBm]");
title("Range estimation");
grid

T1R1Ex = [];
T1R2Ex = [];
T1R3Ex = [];
T2R1Ex = [];
T2R2Ex = [];
T2R3Ex = [];
T3R1Ex = [];
T3R2Ex = [];
T3R3Ex = [];
T4R1Ex = [];
T4R2Ex = [];
T4R3Ex = [];

for(i = 1:32)
    T1R1Ex = horzcat(T1R1Ex,T1R1);
    T1R2Ex = horzcat(T1R2Ex,T1R2);
    T1R3Ex = horzcat(T1R3Ex,T1R3);
    T2R1Ex = horzcat(T2R1Ex,T2R1);
    T2R2Ex = horzcat(T2R2Ex,T2R2);
    T2R3Ex = horzcat(T2R3Ex,T2R3);
    T3R1Ex = horzcat(T3R1Ex,T3R1);
    T3R2Ex = horzcat(T3R2Ex,T3R2);
    T3R3Ex = horzcat(T3R3Ex,T3R3);
    T4R1Ex = horzcat(T4R1Ex,T4R1);
    T4R2Ex = horzcat(T4R2Ex,T4R2);
    T4R3Ex = horzcat(T4R3Ex,T4R3);
end
%% 

T(:,:,1) = DopplerRange(T1R1Ex,R_max,1);
T(:,:,2) = DopplerRange(T1R2Ex,R_max,0);
T(:,:,3) = DopplerRange(T1R3Ex,R_max,0);

T(:,:,4) = DopplerRange(T2R1Ex,R_max,0);
T(:,:,5) = DopplerRange(T2R2Ex,R_max,0);
T(:,:,6) = DopplerRange(T2R3Ex,R_max,0);

T(:,:,7) = DopplerRange(T3R1Ex,R_max,0);
T(:,:,8) = DopplerRange(T3R2Ex,R_max,0);
T(:,:,9) = DopplerRange(T3R3Ex,R_max,0);

T(:,:,10) = DopplerRange(T4R1Ex,R_max,0);
T(:,:,11) = DopplerRange(T4R2Ex,R_max,0);
T(:,:,12) = DopplerRange(T4R3Ex,R_max,0);


%% 

t1 = T(:,:,1);
[row, col] = find(ismember(t1, max(t1(:))));

[~,~,widths,proms] = findpeaks(abs(T(:,col,1)));
thr1 = (sum(maxk(proms,15)))/15;
thr2 = max(widths)*2;
[pks,locs] = findpeaks(abs(T(:,col,1)),'MinPeakProminence',thr1, ...
    'MinPeakDistance',thr2);

for i = 1:length(locs)
    t1=[T(locs(i),col,1) T(locs(i),col,2) T(locs(i),col,3) ...
        T(locs(i),col,4) T(locs(i),col,5) T(locs(i),col,6) ...
        T(locs(i),col,7) T(locs(i),col,8) T(locs(i),col,9) ... 
        T(locs(i),col,10) T(locs(i),col,11) T(locs(i),col,12)];
    
    figure

    Na = length(t1)*10;
    angleFFT = 10*log10(abs(fftshift(fft(t1,Na))));

    rA = linspace(-pi/2,pi/2,Na);
%     [~,~,widths,proms] = findpeaks(angleFFT);
%     thr1 = (sum(mink(proms,1)))*100;
%     thr2 = max(widths)*0.1;
%     [~,locsAng] = findpeaks(angleFFT,'MinPeakProminence',thr1);
    [~,locsAng] = find(ismember(angleFFT, maxk(angleFFT(:),2)));
    anglesRad = rA(locsAng);
    anglesDeg = [];
    for j = anglesRad
       anglesDeg = cat(2,anglesDeg,rad2deg(asin(j/pi)));
    end

    hold on
    plot(rA,angleFFT)
    scatter(rA(locsAng),angleFFT(locsAng));
    title("Spatial spectrums from the virtual array");
    legend("Power [dB]");
    xlabel("Angle [rad]");
    ylabel("Power [dB]");
    xlim([-pi/2,pi/2]);
    grid
    hold off
end
