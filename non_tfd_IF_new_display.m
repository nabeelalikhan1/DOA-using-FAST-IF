

function [fidexmult,TFC] = non_tfd_IF_new_display(Sig,win_length, num, delta,L,thr,Thr)
% SET OF FRACTIONAL WINDOWS
w=gausswin(win_length,1);
l=0:length(Sig)-1;

for iii=1:length(Sig)
    WW(iii,:)=exp(-1i*(iii)*2*pi*l/length(Sig));
end
e_max=0;
%tic;
i=0;
window_rot=zeros(2*L+1,win_length);
for k=-L+1:1:L-1
    i=i+1;
    window_rot(i,:)=fracft(w,0.99* k/L);%0.05
end
%save('window_rot','window_rot');
%load('window_rot');
%toc;
TFC=zeros(length(Sig),length(Sig),2*L+1);
    Sig_extended=[zeros(1,floor(win_length/2)) Sig zeros(1,floor(win_length/2))];

    w_signal=zeros(1,length(Sig));
for iii=1:num
         Siga=filter(ones(1,win_length),1,abs(Sig));
    [~,t_start]=max(Siga(floor(win_length/2)+1:end-floor(win_length/2)));
    t_start=t_start(1)+floor(win_length/2);
    v_m=0;
    for i=1:2*L+1
        FF=abs(fft(Sig(t_start-floor(win_length/2):t_start+floor(win_length/2)).*window_rot(i,:),length(Sig)));
        
        [v,index]=max(FF(1:end/2));
        if v_m<v
            freq_start=index;  %peak_frequency location
            frac_wind_index_start=i;   %chirp rate location
            v_m=v;
        end
    end
    v_oldd=v_m;
    
    
    IF=zeros(1,length(Sig))-1;
    IF(t_start)=freq_start;
    f0=freq_start;
    frac_wind_index=frac_wind_index_start;

    TFC(t_start,freq_start,frac_wind_index_start)=iii;
    
    
    
    for t_index=t_start+1:length(Sig)
        v_m=0;
        k=f0-delta:1:f0+delta;
        
        k(k>length(Sig)/2-1)=length(Sig)/2-1;
        k(k<=0)=1;
        if frac_wind_index<2
            frac_wind_index=2;
        elseif frac_wind_index>2*L-1
            
            frac_wind_index=2*L-1;
        end
        for i=frac_wind_index-1:frac_wind_index+1    % FOR ALL WINDOWS
            w_signal(1:win_length)=(Sig_extended(t_index:t_index+win_length-1).*window_rot(i,:)); %WINDOWED SIGNAL
            for jjj=1:length(k)%jj
                V(jjj)=(abs(sum(w_signal.*WW(k(jjj),:))));
            end
            [v,index]=max(V);
            if v_m<v
                f0=k(index);  %peak_frequency location
                frac_wind_index=i;   %chirp rate location
                v_m=v;
            end
        end
        if v_m<Thr*v_oldd
            
            break;
        end
        % v_old=v_m;
        TFC(t_index,f0,frac_wind_index)=iii;
        %  f0
        IF(t_index)=f0;
        %        IF(t_index)
    end
    
    
    f0=freq_start;
    frac_wind_index=frac_wind_index_start;
    for t_index=t_start-1:-1:1
        v_m=0;
        
        k=f0-delta:1:f0+delta;
        k(k>length(Sig)/2-1)=length(Sig)/2-1;
        k(k<=0)=1;
        if frac_wind_index<2
            frac_wind_index=2;
        elseif frac_wind_index>2*L-1
            
            frac_wind_index=2*L-1;
        end
        
        for i=frac_wind_index-1:frac_wind_index+1    % FOR ALL WINDOWS
            
            
            w_signal(1:win_length)=(Sig_extended(t_index:t_index+win_length-1).*window_rot(i,:)); %WINDOWED SIGNAL
            for jjj=1:length(k)%jj
                V(jjj)=abs(sum(w_signal.*WW(k(jjj),:)));
            end
            [v,index]=max(V);
            if v_m<v
                f0=k(index);  %peak_frequency location
                frac_wind_index=i;   %chirp rate location
                v_m=v;
            end
            
        end
        
        if v_m<Thr*v_oldd
            
            break;
        end
        %v_old=v_m;
        TFC(t_index,f0,frac_wind_index)=iii;
        
        IF(t_index)=f0;
    end
    
    IF=IF/(length(Sig));
    %figure; plot(IF)
    Phase=2*pi*filter(1,[1 -1],IF);
    s_dechirp=exp(-1i*Phase);
    
    
    LL=delta;
    %TF filtering for each sensor
    s1 = Sig.*(s_dechirp);
    s2=fftshift(fft(s1));
    %figure; plot(abs(s2));
    if sum(abs(s2(length(Sig)/2-LL-1:length(Sig)/2+LL-1).^2))<e_max*thr
        break;
    else
        %  e_max
        if sum(abs(s2(length(Sig)/2-LL-1:length(Sig)/2+LL-1).^2))>e_max
            e_max=sum(abs(s2(length(Sig)/2-LL-1:length(Sig)/2+LL-1).^2));
        end
        % Energy of the last component
        s2(length(Sig)/2-LL:length(Sig)/2+LL)=0;
        s2=ifft(ifftshift(s2)).*conj(s_dechirp);
        Sig=s2;%-extr_Sig(iii);
        fidexmult(iii,:) = IF;
    end
    
end
end


