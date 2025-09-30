clear
close all
clc
load ThicknessPF5;
data=E;
[m,n,q] = size(data);
for i=1:m
    for j=1:n
        tmax(i,j)=find(data(i,j,:)==max(data(i,j,:)),1,'first');
        if(max(data(i,j,:))<=0.015)
            tmax(i,j)=0;
        end
    end
end

for i=1:m
    for j=1:n
        if(tmax(i,j)>=1500)
            tmax(i,j)=0;
        end
    end
end

D = tmax.*0.033e-12.*3e8./(1.8-1);

D1 = D(:,35:end);

figure
surf(D)
shading interp