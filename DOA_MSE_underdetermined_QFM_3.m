close all;
clear all;
N_sensors=2;
n=0:127;

%addpath('D:\D\win64_bin\win64_bin');
addpath('E:\tfsa_5-5\windows\win64_bin');


%crossing componentsi8

s1=exp(2*pi*1i*(0.15*n+0.35*n.^3/(128*128*3)));
s2=exp(2*pi*1i*(0.45*n-0.35*n.^3/(128*128*3)));
s3=exp(2*pi*1i*(0.025*n+1*0.35*n.^3/(128*128*3)));

perc=0.4;

t_IF(1,:)=0.05+0*0.35*n/128;
t_IF(2,:)=0.15+0.35*n/128;
t_IF(3,:)=0.45-0.35*n/128;

s = [(s1.') (s2.') (s3.')];

n_sources=3;
N_C=3;
s_orig=s;

% set mixing matrix A
       theta = [-5,0,5]*pi/180;   % sensor separation angles in radians
LL=100;
index=0;
delta=2;
for SNR=-10:5:10
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
        ss= multi_sensor_source_separation_ridge_tracking_m(X, n_sources, delta,N_sensors);
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
        mmssee_ridge_detection_tracking(ii)=mean((sort(y1/10)-sort(theta9/10)).^2);

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
        [ss,IF_out] = multi_sensor_source_separation_spatial_TF_direction(X, n_sources, 3,N_sensors);
          %      [IFF,ss] = Multi_Sensor_FAST_IF(X,N_sensors,65, n_sources, 2,100,0,0);

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
        theta1=theta9*180/pi;
        y1=y1-90;
        mmssee_sadtfd(ii)=mean((sort(y1/10)-sort(theta9/10)).^2);
%         
        
        %tic
        D   = mtfd(X, 'ckd',1, 0.3, 0.3, length(X));
        %toc
        %%% Averaged Auto-TFD
        D_avg = zeros(length(X), length(X));
        for mm = 1:N_sensors, D_avg = D{mm,mm} + D_avg; end
        D_avg = D_avg./N_sensors;
        %%% Selection of high-energy (t,f) points
        thr = perc*max(max(D_avg));
        Tr = abs(D_avg) >= thr;
        [F_trace, ~] = find(Tr);
        n_p = length(F_trace);
        D_s = zeros(N_sensors, N_sensors);
        for m1 = 1:N_sensors
            for m2 = 1:N_sensors
                D_s(m1,m2) = (1/n_p).*sum(sum(D{m1,m2}.*Tr));
            end
        end
        theta1=-90:1:90;

        %%% DOA Estimation
        P_tf_music_ckd = tf_music(D_s, n_sources, N_sensors, 2,1, theta1);
       % P_tf_music_ckd =TMMUSIC(D_s, 2, N_sensors, n_sources, 1, theta1);
        [~,y1]=findpeaks(P_tf_music_ckd,'NPeaks',4,'MinPeakDistance',10,'Threshold',0.001);
        y1=y1-90;
        %figure;plot(P_tf_music_ckd)
        if length(y1)<n_sources
            y1(length(y1)+1:n_sources)=0;
        elseif length(y1)>n_sources
            y1=y1(1:n_sources);
        end
        
                
        
        mmssee_tf_music(ii)=mean((sort(y1/10)-sort(theta9/10)).^2);

    end
    index=index+1;
    %mean(mmssee)
    snr_mse_fastIF(index)=mean(mmssee_FASTIF)
    snr_mse_sadtfd(index)=mean(mmssee_sadtfd)
    snr_mse_tf_music(index)=mean(mmssee_tf_music)
    snr_mse_ridge_tracking(index)=mean(mmssee_ridge_detection_tracking)
end

%SNR=-10:2:0;
SNR=-10:5:10;
% FROM SIMULATION DONE in "Novel direction of arrival estimation using
% spatial adaptive"
%snr_mse_sadtfd=[0.4418    0.0994    0.0615    0.0414    0.0272    0.0204];
plot(SNR,10*(log10(snr_mse_sadtfd)),'--md','linewidth',2);
hold on;
plot(SNR,10*(log10(snr_mse_fastIF)),'r','linewidth',3);
hold on;
plot(SNR,10*(log10(snr_mse_tf_music)),'g','linewidth',3);
hold on;
%plot(SNR,10*(log10(snr_mse_post_proc)),'y','linewidth',2);
%hold on;
plot(SNR,10*(log10(snr_mse_ridge_tracking)),'b:','linewidth',3);

%openfig('DOA_MSE_over_determined -10 db to 10 db.fig')
xlabel('Signal to Noise Ratio');
ylabel('Mean Square Error (dB)');
legend('Adaptive Spatial TFDs','The FAST IF','Time-frequency Music','DOA based on IF estimation using ridge tracking');
%legend('The Proposed Method','Time-frequency Music','DOA based on IF estimation using ridge tracking');
