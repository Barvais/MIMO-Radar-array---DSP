clc
clear
close all

%Define constants and importing raw data.
c = physconst('lightspeed');
BW = 2e9;

Mag=load('S12Mag.tab');
Phase=load('S12Phase.tab');

%Extracting relevant data for processing.
S21=Mag(:,2).*exp(j*Phase(:,2));


N=length(S21);
Delta=c/(2*BW);
R_max=(N-1)*Delta;

%Performing inverse FFT to extract the range from the S(TX,RX).
X = ifft(S21,N*10);

r=linspace(0,R_max,N*10);

%Finding peaks while ignoring small peaks while using avarage on the top 15
%largest prominences.
[~,~,widths,proms] = findpeaks(abs(X));
thr1 = (sum(maxk(proms,15)))/15;
thr2 = max(widths)*2;
[pks,locs] = findpeaks(abs(X),'MinPeakProminence',thr1,'MinPeakDistance',thr2);

%Printing range plot in dB.
figure;
hold on;
plot(r,10*log10(abs(X)));
X1 = 10*log10(abs(X));
scatter(r(locs),X1(locs));
hold off;
xlabel("Range [m]");
ylabel("Power [dBm]");
title("Range estimation TX1-RX1 ");
legend("Power [dBm]");
xlim([0 200]);
grid

