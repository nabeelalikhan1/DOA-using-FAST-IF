clc
clear
close all
mul=1;
SampFreq = 128*mul;
addpath('D:\tfsa_5-5\windows\win64_bin');
t = 0:1/SampFreq:1-1/SampFreq;
Sig1 = 1*exp(mul*1i*(1*pi*(30*t.^3))+mul*1i*(2*pi*(0*t))); %300t或者150t
Sig2 = 1*exp(mul*1i*(-1*pi*(30*t.^3))+mul*1i*(1*pi*(90*t))); %300t或者150t

Sig3 = exp(mul*1i*(1*pi*(30*t +30*t.^3)));
Sig4 =1*exp(mul*1i*(1*pi*(120*t -30*t.^3)));
num=2;
Sig =1*Sig1.*1 +	0*Sig4 +0*Sig3+1*Sig2.*1;
Sig(2,:)=Sig;
Sig(3,:)=Sig(1,:);

Sig=awgn(Sig,-5,'measured');
IF_O(:,1)=mul*90*t.^2/2;
IF_O(:,2)=-mul*90*t.^2/2+mul*90/2;



%[fidexmult,A] = non_tfd_IF_new_display(Sig,length(Sig)/(2)-1, num, 2,100,0,0);

[fidexmult,Xout,A] = Multi_Sensor_FAST_IF(Sig,3,length(Sig)/(2)-1, num, 2,100,0,0);
figure;
plot(t,IF_O,':',t,SampFreq*fidexmult,'--','linewidth',3);
axis([0 1 0 64]);
xlabel('Time (s)');
ylabel('Instantaneous Frequency (Hz)')
figure;idx = find(A>=1);
[X, Y, Z] = ind2sub(size(A), idx);
pointsize = 30;
%scatter3(X(:), Y(:), Z(:), pointsize, A(idx));
scatter3(Z(:), X(:), Y(:), pointsize, A(idx));
xlabel('Frequency (Hz)');
ylabel('Time(s)');
zlabel('Chirp rate ');

A=sum(A,3);
figure;
tfsapl(Sig,A)

% HADTFD BASED
