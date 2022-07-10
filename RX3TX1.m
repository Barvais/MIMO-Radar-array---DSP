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

N=length(TX1RX123)

Delta=c/(2*BW);
R_max=(N-1)*Delta

r=linspace(0,R_max,N*10);
rA=linspace(-pi/2,pi/2,N*10);

T1R1 = ifft(TX1RX1,N*10);
T1R2 = ifft(TX1RX2,N*10);
T1R3 = ifft(TX1RX3,N*10);

T1R1Ex = [];
T1R2Ex = [];
T1R3Ex = [];
for(i = 1:128)
    T1R1Ex = horzcat(T1R1Ex,T1R1);
    T1R2Ex = horzcat(T1R2Ex,T1R2);
    T1R3Ex = horzcat(T1R3Ex,T1R3);
end

% [~,~,widths,proms] = findpeaks(abs(T1R1));
% thr1 = (sum(maxk(proms,15)))/15;
% thr2 = max(widths)*2;
% [pks,locs] = findpeaks(abs(T1R1),'MinPeakProminence',thr1,'MinPeakDistance',thr2);

T1(:,:,1) = T1R1Ex;
T1(:,:,2) = T1R2Ex;
T1(:,:,3) = T1R3Ex;
%% 

DR(:,:,1)=DoplerRange(T1R1Ex,R_max);
DR(:,:,2)=DoplerRange(T1R2Ex,R_max);
DR(:,:,3)=DoplerRange(T1R3Ex,R_max);
d = DR(:,:,1);
[row, col] = find(ismember(d, max(d(:))));
% row = peaks(tAbs(:,col))

t1=[DR(row,col,1) DR(row,col,2) DR(row,col,3)];
Na = length(t1);
% test = [t1 t1 t1 t1 t1 t1 t1 t1 t1 t1];
rA = linspace(-pi/2,pi/2,Na);
plot(rA,abs(fftshift(fft(t1,Na))))

%% 
t = fftn(T1);
[~,I] = max(t);
% I=locs;
% t2=[t(I(2),:,1) t(I(2),:,2) t(I(2),:,3)]; 
t1=[t(I(1),:,1) t(I(1),:,2) t(I(1),:,3)]; 
%     t(I(2),:,4) t(I(2),:,5) t(I(2),:,6) ...
%     t(I(2),:,7) t(I(2),:,8) t(I(2),:,9) ... 
%     t(I(2),:,10) t(I(2),:,11) t(I(2),:,12)
% t2 = angle(t1);
close all
figure
% test = [t1 t1 t1 t1 t1];
Na = length(t1);
rA = linspace(-pi/2,pi/2,Na);
plot(rA,abs((fft(t1,Na))))
% figure
% plot(rA,abs((fft(t2,Na))))

a1=rad2deg(asin((angle(T1R1(I(1)))-angle(T1R2(I(1))))/pi))
a2=rad2deg(asin((angle(T1R2(I(1)))-angle(T1R3(I(1))))/pi))
a3=rad2deg(asin((angle(T1R1(I(1)))-angle(T1R3(I(1))))/(2*pi)))
a4=rad2deg(asin((angle(T1R1(I(2)))-angle(T1R2(I(2))))/pi))
a5=rad2deg(asin((angle(T1R2(I(2)))-angle(T1R3(I(2))))/pi))
a6=rad2deg(asin((angle(T1R1(I(2)))-angle(T1R3(I(2))))/(2*pi)))

(a1+a2+a3)/3
(a4+a5+a6)/3
%% 

% test = T1(locs,:,:)
% test1 = [test(:,:,1) test(:,:,2) test(:,:,3)]
%% 

t = fftn(T1);
t1 = angle(t);
abst = abs(t);
a = [abst(:,:,1) abst(:,:,2) abst(:,:,3)];
a1 = [t1(:,:,1) t1(:,:,2) t1(:,:,3)];

plot(rA,abst(:,:,1))
%% 

% figure
% imagesc(10*log10(a));

figure;
subplot(2,1,1);
plot(r,10*log10(abs(T1R1)));
xlabel("Range [m]");
ylabel("P [dBm]");
title("Range estimation TX1 to RX1");
grid
subplot(2,1,2);
plot(r,angle(T1R1));

figure;

plot(r,10*log10(abs(T1R2)));
xlabel("Range [m]");
ylabel("P [dBm]");
title("Range estimation TX1 to RX2");
grid

figure;

plot(r,10*log10(abs(T1R3)));
xlabel("Range [m]");
ylabel("P [dBm]");
title("Range estimation TX1 to RX3");
grid

% rad2deg(asin((-1.6099+1.5897)/pi))

% X1 = 10*log10(abs(X));
%scatter(r(locs),X1(locs));

% xlabel("Range [m]");
% ylabel("P [dBm]");
% title("Range estimation");
% grid