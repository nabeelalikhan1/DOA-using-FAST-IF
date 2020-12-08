close all;
clear all;
N_sensors=3;
n=0:127;

%addpath('D:\D\win64_bin\win64_bin');
%addpath('D:\tfsa_5-5\windows\win64_bin');

%crossing componentsi8

s1=exp(2*pi*1i*(0.05*n+0.45*n.^3/(128*128*3)));
s2=exp(2*pi*1i*(0.125*n+0.2*n.^3/(128*128*3)));
s2=exp(2*pi*1i*(0.45*n-0.45*n.^3/(128*128*3)));

perc=0.4;

t_IF(1,:)=0.05+0.45*3*n.^2/(128*128*3);
t_IF(2,:)=0.125+0.2*3*n.^2/(128*128*3);
t_IF(2,:)=0.45-0.45*3*n.^2/(128*128*3);

%s = [(s1.') (s2.') (s3.')];
s = [(s1.') (s2.') ];

n_sources=2;
N_C=2;
s_orig=s;

% set mixing matrix A
%theta = [15,30,50]*pi/180;   % sensor separation angles in radians
        theta = [-12,12]*pi/180;   % sensor separation angles in radians
%        theta = [-12,0,12]*pi/180;   % sensor separation angles in radians
       theta = [-5,5]*pi/180;   % sensor separation angles in radians
LL=100;
index=0;
delta=2;
for SNR=-10:2:10
 %for SNR=-2:2:0  
    for ii=1:LL
        
        
        A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));  % mixing matrix A
        
        
        X = A*s.';                             % mixed source
        theta9=round(theta *180/pi);
        % generate noise
        
        sigma = 10^(-SNR/20);
        w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise
        
        X=X+w;
        
      %  [ss] = multi_sensor_source_separation_using_TRUE_IF(t_IF,X,N_C, 3,N_sensors);
      %tic
       
        %tic
        [IFF,ss] = Multi_Sensor_FAST_IF(X,N_sensors,65, n_sources, 2,100,0,0);
        
        %toc
        
        for iii=1:n_sources
            for jjj=1:N_sensors
                a(jjj,:)=ss(jjj,iii,:);
            end
            theta1=-90:1:90;
            
            p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta1');
            [x,y]=max(p);
            y1(iii)=y;
        end
        
        y1=y1-90;
        
        mmssee_FASTIF(ii)=mean((sort(y1/10)-sort(theta9/10)).^2);
        
        
        %original code commented
        %[ss,IF_out] = multi_sensor_source_separation_spatial_TF_direction(X, n_sources, 3,N_sensors);
       
    end
    index=index+1;
    %mean(mmssee)
    snr_mse_fastIF(index)=mean(mmssee_FASTIF)
    
end

%SNR=-10:2:0;
SNR=-10:2:10;
% FROM SIMULATION DONE in "Novel direction of arrival estimation using
% spatial adaptive"
%snr_mse_sadtfd=[0.4418    0.0994    0.0615    0.0414    0.0272    0.0204];

openfig('DOA_MSE_over_determined -10 db to 10 db.fig');
hold on;

plot(SNR,10*(log10(snr_mse_fastIF)),'r','linewidth',3);

%xlabel('Signal to Noise Ratio');
%ylabel('Mean Square Error (dB)');
%legend('Adaptive Spatial TFDs','The FAST IF','Time-frequency Music','DOA based on IF estimation using ridge tracking');
%legend('The Proposed Method','Time-frequency Music','DOA based on IF estimation using ridge tracking');
