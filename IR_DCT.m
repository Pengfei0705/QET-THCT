clear
close all
clc
load ThicknessPF5;
data=ImageLabHSI.dataCube;
[a,b,c]=size(data);
t=ImageLabHSI.axisCoordsL;
fr=1/0.033;    % frame rate
f = linspace(0,fr,c);
sf = f(c/4:-1:2);
sf(1:499)=[];
t=linspace(0,1/fr*c,c);
for i=1:length(sf)
   ref(i,:) = cos(2*pi*sf(i)*2*t);
   quad(i,:) = sin(2*pi*sf(i)*2*t);
end

figure
plot(t,ref(end,:))
hold on
plot(t,ref(249,:))
plot(t,quad(end,:))
plot(t,quad(249,:))

export = [t;ref(249:250,:);quad(249:250,:)]';






for i=1:length(sf)
    FInp(i,:) = conj(fft(ref(i,:)));
    FQuad(i,:) = conj(fft(quad(i,:)));
end

s=reshape(data(104,79,:),1,c);
for k=1:length(sf)
S0=ifft(FInp(k,:).*fft(s));
S1=ifft(FQuad(k,:).*fft(s));
SS0(k,:) = S0;
SS1(k,:) = S1;
Amp(k,:) = sqrt(S0.^2+S1.^2);
Pha(k,:) = atan2(S1,S0); 
end

figure
plot(SS0(1,:))
hold on
plot(SS1(1,:))


figure
plot(Amp(1,:))
hold on
plot(Amp(100,:))
plot(Amp(200,:))

for k=1:length(sf)
    ref(k,:) = rescale(ref(k,:),0,1);
    Amp(k,:) = rescale(Amp(k,:),0,1);
    Pha(k,:) = rescale(Pha(k,:),0,1);
end

figure
Y = SS0(170,:)';
plot(Y)

figure
Y = ref(183,:)';
plot(Y)

figure
imagesc(ref)

tt= 0:0.1:100;
figure
y= 100*sin(tt)./tt;
plot(tt,y)

figure
imagesc(Pha)

xlabel('time (s)')
ylabel('CC Amplitude')
Amp=Amp';



%%
clear
close all
clc

data=ImageLabHSI.dataCube;

[a,b,c]=size(data);


fr=1/0.033;    % frame rate
f = linspace(0,fr,c);
sf = f(c/4:-1:2);
sf(1:599)=[];
t=linspace(0,1/fr*c,c);
for i=1:length(sf)
   ref(i,:) = cos(2*pi*sf(i)*t);
   quad(i,:) = sin(2*pi*sf(i)*t);
end

for i=1:length(sf)
    FInp(i,:) = conj(fft(ref(i,:)));
    FQuad(i,:) = conj(fft(quad(i,:)));
end

se = floor(c/2);
temp_Amplitude1 = zeros(a, b, length(sf));
temp_Phase1 = zeros(a, b, length(sf));

parfor m = 1:a
    for n = 1:b
        temp_Amp = zeros(1, length(sf));
        temp_Pha = zeros(1, length(sf));
        
        for i = 1:length(sf)
            s = reshape(data(m, n, :), 1, c);
            S0 = ifft(FInp(i, :) .* fft(s));
            S1 = ifft(FQuad(i, :) .* fft(s));
            Amp = sqrt(S0.^2 + S1.^2);
            Pha = atan2(S1, S0);
            temp_Amp(i) = max(Amp(1:se));
            k = find(Amp == max(Amp(1:se)), 1, 'first');
            temp_Pha(i) = Pha(k);
        end
        
        temp_Amplitude1(m, n, :) = temp_Amp;
        temp_Phase1(m, n, :) = temp_Pha;
    end
end

% Combine temporary arrays into final result
Amplitude1 = temp_Amplitude1;
Phase1 = temp_Phase1;

yy1 = fliplr(reshape(Amplitude1(60,60,:),[],150));
yy2 = reshape(FFT_Amp_1(60,60,:),[],150);

figure
plot(yy1,'r-')
hold on
plot(yy2,'k-')


