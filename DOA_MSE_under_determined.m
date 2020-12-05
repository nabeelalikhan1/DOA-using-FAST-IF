close all;
clear all;
N_sensors=2;
n=0:127;

%addpath('D:\D\win64_bin\win64_bin');
addpath('D:\tfsa_5-5\windows\win64_bin');
%addpath('E:\Published Papers\DOA estimation of intersecting components 2018\Matlab code');
%addpath('E:\Published Papers\DOA ESTIMATION VITERBI\Multi-sensor IF estimation code');

%crossing componentsi8

s1=exp(2*pi*1i*(0.05*n+0.45*n.^3/(128*128*3)));
s2=exp(2*pi*1i*(0.45*n-0.45*n.^3/(128*128*3)));

perc=0.4;



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
LL=200/2;
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
        clear a;
        clear ss;
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
        
        
    %    [ss,IFF] = multi_sensor_source_separation_spatial_diversity_post_proc(X, N_C, 3,N_sensors);
        %figure; plot(IFF');
        
     %   for iii=1:n_sources
      %      for jjj=1:N_sensors
       %         a(jjj,:)=ss(jjj,iii,:);
        %    end
         %   theta1=-90:1:90;
            
          %  p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta1');
           % [x,y]=max(p);
           % y1(iii)=y;
        %end
        
        %y1=y1-90;
        
        %mmssee_post_proc(ii)=mean((sort(y1/10)-sort(theta9/10)).^2);
        %tic
        %original code commented
        [ss,IF_out] = multi_sensor_source_separation_spatial_TF_direction(X, n_sources, 3,N_sensors);
       % [IFF,ss] = Multi_Sensor_FAST_IF(X,N_sensors,65, n_sources, 2,100,0,0);

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
        

    end
    index=index+1;
    %mean(mmssee)
    snr_mse_fastIF(index)=mean(mmssee_FASTIF)
    snr_mse_sadtfd(index)=mean(mmssee_sadtfd)
    snr_mse_ridge_tracking(index)=mean(mmssee_ridge_detection_tracking)
end

%SNR=-10:2:0;
SNR=-10:2:10;
%SNR=0:2:10;
% FROM SIMULATION DONE in "Novel direction of arrival estimation using
% spatial adaptive"
plot(SNR,10*(log10(snr_mse_sadtfd)),'--md','linewidth',2);
hold on;
plot(SNR,10*(log10(snr_mse_fastIF)),'r','linewidth',3);
hold on;

%plot(SNR,10*(log10(snr_mse_post_proc)),'y','linewidth',2);
%hold on;
plot(SNR,10*(log10(snr_mse_ridge_tracking)),'b:','linewidth',3);

xlabel('Signal to Noise Ratio');
ylabel('Mean Square Error (dB)');
legend('SADTFD-Viterbi based DOA estimation', 'The Proposed Method','DOA based on IF estimation using ridge tracking');

%legend('The Proposed Method','DOA based on IF estimation using ridge tracking');
