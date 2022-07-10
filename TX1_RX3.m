clc
clear
close all

c = physconst('lightspeed');
BW = 2e9;
MagTX1=load('MagsTX1.tab');
PhaseTX1=load('PhasesTX1.tab');

TX1RX123=MagTX1(:,2:end).*exp(j*PhaseTX1(:,2:end));


TX1RX1=TX1RX123(:,1);
TX1RX2=TX1RX123(:,2);
TX1RX3=TX1RX123(:,3);

N=length(TX1RX123);

Delta=c/(2*BW);
R_max=(N-1)*Delta;

r=linspace(0,R_max,N*10);
rA=linspace(-pi/2,pi/2,12);

T1R1 = ifft(TX1RX1,N*10);
T1R2 = ifft(TX1RX2,N*10);
T1R3 = ifft(TX1RX3,N*10);

figure

plot(r,10*log10(abs(T1R1)));
xlabel("Range [m]");
ylabel("P [dBm]");
title("Range estimation");
grid

T1R1Ex = [];
T1R2Ex = [];
T1R3Ex = [];

for(i = 1:32)
    T1R1Ex = horzcat(T1R1Ex,T1R1);
    T1R2Ex = horzcat(T1R2Ex,T1R2);
    T1R3Ex = horzcat(T1R3Ex,T1R3);
end

%% 

T(:,:,1) = DopplerRange(T1R1Ex,R_max,1);
T(:,:,2) = DopplerRange(T1R2Ex,R_max,0);
T(:,:,3) = DopplerRange(T1R3Ex,R_max,0);


%% 

t1 = T(:,:,1);
[row, col] = find(ismember(t1, max(t1(:))));

[~,~,widths,proms] = findpeaks(abs(T(:,col,1)));
thr1 = (sum(maxk(proms,15)))/15;
thr2 = max(widths)*2;
[pks,locs] = findpeaks(abs(T(:,col,1)),'MinPeakProminence',thr1, ...
    'MinPeakDistance',thr2);

for i = 1:length(locs)
    t1=[T(locs(i),col,1) T(locs(i),col,2) T(locs(i),col,3)] ;
    
    figure

    Na = length(t1)*10;
    angleFFT = 10*log10(abs(fftshift(fft(t1,Na))));

    rA = linspace(-pi/2,pi/2,Na);
    [~,locsAng] = findpeaks(angleFFT);
    
    anglesRad = rA(locsAng);
    anglesDeg = [];
    for j = anglesRad
       anglesDeg = cat(2,anglesDeg,rad2deg(asin(j/pi)));
    end
    
    hold on
    plot(rA,angleFFT)
    scatter(rA(locsAng),angleFFT(locsAng));
    title("Spatial spectrums from the physical array");
    legend("Power [dB]");
    xlabel("Angle [rad]");
    ylabel("Power [dB]");
    xlim([-pi/2,pi/2]);
    grid
    hold off
end


