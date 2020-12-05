function [Inew,D,orient]= SADTFD(S,alpha1,alpha2,N)

%length

%alpha=0.25;
R=3;
D = mtfd(S, 'WVD');%,0.04,0.04);

[MM,NN]=size(D);
D_avg = zeros(length(S), length(S));
for mm = 1:MM, D_avg = D{mm,mm} + D_avg; end
M=180;
B=ones(5,5);
D_avg=filter2(B,D_avg);


[X,Y]=meshgrid(-1:2/N:1,-1:2/N:1);

for kk=0:M/R-1
    angle=pi*kk*R/M;
    
    X1=X*cos(angle)-Y*sin(angle);
    Y1=X*sin(angle)+Y*cos(angle);
    
    
    A=exp((-1/2)*(((alpha1*X1).^2)+(alpha2*Y1).^2));
    %  BB1(kk+1,:,:)=exp((-1/2)*(((alpha1*X1).^2)+(2*alpha2*Y1).^2));
    A=A.*(1-alpha2*alpha2*Y1.^2);
    
    A=A/sum(sum(abs(A)));
    BB{kk+1} = A;
    
end

PQ = paddedsize(size(D_avg));
F = fft2(double(D_avg), PQ(1), PQ(2));

FIabs = fft2(double(abs(D_avg)), PQ(1), PQ(2));

jjjj=0;
for ii=0:M/R-1
    %     B(:,:)=BB(ii+1,:,:);
    %  if or(ii*R<70,ii*R>110)
    %        ii
    %       BB{1}
    B = BB{ii+1};
    
    %B1(:,:)=BB1(ii+1,:,:);
    H = fft2(double(B), PQ(1), PQ(2));
    H1 = fft2(double(B), PQ(1), PQ(2));
    
    F_fH = H.*F;
    ffi = real(ifft2(F_fH));
    II3(:,:,jjjj+1)=(ffi(round(length(B)/2):end-round(length(B)/2), round(length(B)/2):end-round(length(B)/2)));
    %F_fH = H.*FIabs;
    
    F_fH = H1.*FIabs;
    
    ffi =( abs(ifft2(F_fH))).^2;
    
    II(:,:,jjjj+1)=(ffi(round(length(B)/2):end-round(length(B)/2), round(length(B)/2):end-round(length(B)/2)));
    %  end
    jjjj=jjjj+1;
end

Inew=zeros(size(ffi(1:end/2, 1:end/2)));
[M1,N1]=size(ffi(1:end/2, 1:end/2));
[b,a]=max(II,[],3);
lag=0;
for m=1:M1+lag
    for n=1:N1+lag
        %a=find(max(II(m,n,:))==II(m,n,:));
        
        %aa=find(min(II(m,n,:))==II(m,n,:));
        %orient(m,n)=aa(1);  % Direction of each pixel
        Inew(m,n)=II3(m,n,a(m,n));
        orient(m,n)=a(m,n);
    end
end




for i=1:MM
    for j=i:NN
        I=D{i,j};
        B=ones(5,5);

        I=filter2(B,I);

        F = fft2(double(I), PQ(1), PQ(2));

        jjjj=0;
        for ii=0:M/R-1
            B = BB{ii+1};
            
            H = fft2(double(B), PQ(1), PQ(2));
            F_fH = H.*F;
            ffi = (ifft2(F_fH));
            II3(:,:,jjjj+1)=(ffi(round(length(B)/2):end-round(length(B)/2), round(length(B)/2):end-round(length(B)/2)));
            jjjj=jjjj+1;
        end
        for m=1:M1+lag
            for n=1:N1+lag
                IS(m,n)=II3(m,n,a(m,n));
            end
        end
        D{i,j}=IS;  
        if i~=j
            D{j,i}=conj(IS);
        end
    end
end


%end