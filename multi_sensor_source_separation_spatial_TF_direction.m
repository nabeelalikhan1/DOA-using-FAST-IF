function [ss,IF_out] = multi_sensor_source_separation_spatial_TF_direction(X, num, delta,n_sensors)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%  output   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fidexmult: obtained ridge index for multi-component signals
for i = 1:num
 %while 1   
  %   i=i+1;
    %Spec=tfr_stft_high(Sig);
    [ADTFD,D,orienttfd]= SADTFD(X,2,20,84);

orient=orienttfd;
%orient=cell(256,256);
for ii=1:length(X)
    for jj=1:length(X)
        Ds=zeros(n_sensors,n_sensors);
        for mm = 1:n_sensors
            for nn=1:n_sensors
                Ds(mm,nn)=D{mm,nn}(ii,jj);
            end
        end
        P = tf_music(Ds, 1, n_sensors, 2, 1, (-90:1:90)');
        [~,orient(ii,jj)]=max(P);
%        [vv,~]  = svd(Ds);
 %       orient{ii,jj}=vv(:,1)/norm(vv(:,1));
    end
end
orient=orient-90;
     c = findridges_new_viterbi_M_adtfd(ADTFD,orienttfd,orient);

    IF=(c)/(2*length(X));
    
    Phase=2*pi*filter(1,[1 -1],IF);
    s_dechirp=exp(-1i*Phase);
    
    %im_label2=bwmorph(im_label2,'dilate',3);
    
    % For each sensor do the following steps
    IF=(c)/(2*length(X));
    IFF(i,:)=IF;
    L=delta;
    %TF filtering for each sensor
    for iii=1:n_sensors
        
        s=(X(iii,:));
        
        %TF filtering for each sensor
        s1 = s.*(s_dechirp);
        s2=fftshift(fft(s1));
        s3=zeros(1,length(s));
        s3(length(s)/2-L:length(s)/2+L)=s2(length(s)/2-L:length(s)/2+L);
        s2(length(s)/2-L:length(s)/2+L)=0;
        s1=ifft(ifftshift((s3))).*conj(s_dechirp);
        s4=ifft(ifftshift((s2))).*conj(s_dechirp);
        X(iii,:)=s4;
        ss(iii,i,:)=s1;
        % nly for visualization    TFD of each extracted components
        %        l=l+tfrwv(conj(s1'));
    end
    
%     s2=fftshift(fft(s1));
%     s3=zeros(1,length(Sig));
%     s3(128-L:128+L)=s2(128-L:128+L);
%     s2(128-L:128+L)=0;
%     extr_Sig=ifft(ifftshift(s3)).*conj(s_dechirp);
%     s2=ifft(ifftshift(s2)).*conj(s_dechirp);
%     
%     %Sig(iii)=Sig(iii)-extr_Sig(iii);
%     Sig=s2;%-extr_Sig(iii);
aaa=ss(:,i,:).^2;
    E(i)=sum(abs(aaa(:)));
    if E(end)<0.2*max(abs(E))
        
    
        break;
    else
        IF_out=IFF;
        ss_out=ss;
    end
    % extr_Sig1=extr_Sig1+extr_Sig;
    %hold on; quiver(1:256,1:256,cos(orienttfd*pi/180),sin(orienttfd*pi/180));
    %fidexmult(i,:) = c;
end

end