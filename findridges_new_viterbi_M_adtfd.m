function IF = findridges_new_viterbi_M_adtfd(Spec,Specangle,orient)
Path=zeros(size(Spec));
Pathweight=zeros(1,length(Spec));
%Specangle(Specangle>90)=-180+Specangle(Specangle>90);
Path(:,1)=1:length(Path);
for i=2:length(Spec)
    for ii=1:length(Spec) % all members of
        Pathw=zeros(1,length(Spec));
        % For transition from ii to j.
        
        for jj=1:length(Spec)
            rr=abs(jj-ii)-2;
            if rr<0
                rr=0;
                
            end
           % angg1=1/(0.01+(sum(orient{jj,i-1}.*orient{ii,i}))/(1));
            angg1=5*abs(orient(jj,i-1)/(i-1)-orient(ii,i));
    %    Pathw(jj)=rank(Spec(:,i),ii)+Pathweight(jj)+1000*rr+1*abs(angg1-1);% weight of ii +
        Pathw(jj)=rank(Spec(:,i),ii)+Pathweight(jj)+1000*rr+1*abs(angg1);% weight of ii +
        
    end
    [value,index]=min(Pathw);
    Pathweight1(ii)=value;
    Path(ii,i)=index;
    orient(ii,i)=orient(index,i-1)+orient(ii,i);
    %orient{ii,i}=orient{index,i-1}+orient{ii,i};
    %orient{ii,i}=orient{ii,i}/norm(orient{ii,i});
end
Pathweight=Pathweight1;
end

IF(1,length(Spec))=0;
[~,index]=min(Pathweight);
for i=length(Spec):-1:1
    IF(i)=index;
    index=Path(index,i);
end
% Code for decoding
end
function a= rank(A,b)
As=sort(A,'descend');
ind=find(As==A(b));
%if A(b)~=0
a=ind(1)-1;
%else
%   a=128;

%end
end