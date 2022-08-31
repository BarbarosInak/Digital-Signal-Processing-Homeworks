%% 1) Butterworth Filter Design

[B,A]=butter(7,0.29,'low');

%% 2) Kaiser Window Design

delta=10^-5;
N=-20*log10(delta);
beta=0.1102*(N-8.7);

M=(N-8)/(2.285*0.15*pi);

M=round(M);

alfa=M/2;

h=zeros(M,1);

for i=1:M+1;
    
    h(i)=(sin(0.325*pi*((i-1)-alfa))/(pi*((i-1)-alfa)))*(besselj(0,beta*(1-((i-1-alfa)/alfa)^2)^(1/2)))/besselj(0,beta);
    
end



%% 3) FIR Filter Design

firpm_filter=firpm(M,[0 0.25 0.4 1] ,[1 1 0 0]);

%% 4) Plotting Designed Filters

% Scaled magnitude responses

figure;
[H1 W1]=freqz(B,A);
plot(W1/pi,abs(H1));
ylabel('Magnitude');
xlabel('Frequency');
title('Butterworth');

figure;
[H2 W2]=freqz(h,1);
plot(W2/pi,abs(H2));
title('Kaiser Window');

figure;
[H3 W3]=freqz(firpm_filter,1);
plot(W3/pi, abs(H3));
title('FIR');

% Phase Responses

figure;
subplot(3,1,1);
phasez(B,A);
title('Phase Response of Butterworth');

subplot(3,1,2); 
phasez(h);
title('Phase Response of Kaiser Window');

subplot(3,1,3); 
phasez(firpm_filter);
title('Phase Response of FIR');

% Impulse Response

figure;
subplot(3,1,1); 
impz(B,A);
title('Impulse Response of Butterworth'); 
xlabel('Sample'); 
ylabel('Amplitude');

subplot(312); 
impz(h,1);
title('Impulse Response of Kaiser Window'); 
xlabel('Sample'); 
ylabel('Amplitude');

subplot(313); 
impz(firpm_filter,1);
title('Impulse Response of FIR'); 
xlabel('Sample'); 
ylabel('Amplitude');


%% 5) High Pass

[num_1,den_1] = iirlp2hp(B, A, 0.25, 0.70);

figure;
freqz(num_1,den_1);

%% 6) Band Pass

[num_2, den_2] = iirlp2bp(B, A, 0.25, [0.3, 0.5]);

figure;
freqz(num_2,den_2);

tf(num_2,den_2,2000,'variable','z^-1')

%% 7) Filter

delta=10^-5;

Nmb=2;

G=[4*10^-3 2*10^-3.75];

w=[0.4*pi 0.7*pi];

N1=-20*log10(G(1));
beta1=0.1102*(N1-8.7);

M1=(N1-8)/(2.285*0.07*pi);

M1=round(M1);

alfa1=M1/2;




N2=-20*log10(G(2));
beta2=0.1102*(N2-8.7);

M2=(N2-8)/(2.285*0.07*pi);

M2=round(M2);

alfa2=M2/2;

h2=zeros(M2,1);

for i=1:M2+1;
    
    h2(i)=((G(2)-delta)/(G(1)-G(2)))*(sin(w(2)*((i-1)-alfa2))/(pi*((i-1)-alfa2)))*(besselj(0,beta2*(1-((i-1-alfa2)/alfa2)^2)^(1/2)))/besselj(0,beta2);
    
end

h1=zeros(M2+1,1);

for i=1:M1+1;
    
    h1(i)=1*(sin(w(1)*((i-1)-alfa))/(pi*((i-1)-alfa)))*(besselj(0,beta1*(1-((i-1-alfa)/alfa)^2)^(1/2)))/besselj(0,beta1);
    
end


hn=h1+h2;

figure;
[H4 W4]=freqz(hn,1);
H5=[H4' flip(H4)'];
plot(linspace(0,2,length(H5)),abs(H5));
ylabel('Magnitude');
xlabel('w rad/sec (x pi)');
title('Magnitude Response');


%% 8)
n = 1:15;
x = cos(2*pi*0.3*n);
factor = 4;

figure
stem(x);
title("Old Plot");
xlabel('n');
ylabel('Magnitude');



x_new= x;
fact = factor;
padlen = length(x_new)*(fact - 1);
zlen = ceil((length(x_new)+1)/2);
z = fft(x_new);

zp = [z(1:zlen) zeros(1, padlen) z(zlen+1:end)];
if  (mod(length(x_new), 2)==0)
    zp(zlen) = zp(zlen)/2; zp(zlen+padlen) = zp(zlen);
end

x_p = real(ifft(zp))*m; 

figure
stem(x_p);
title("New Plot");
xlabel('n');
ylabel('Magnitude');

%%
teta_p=0.25*pi;
w_p=0.7*pi;
ts=9800;

alfa_n=sin((teta_p-w_p)/2)/sin((teta_p+w_p)/2);

B_p=[-alfa_n 1];
A_p=[1 -alfa_n];

Z_p=tf(B_p,A_p,ts,'variable','z^-1');

ans=0.001836*((1+Z_p)^4)/((1-1.5548*Z_p+0.6493*Z_p^2)*(1-1.4996*Z_p+0.8482*Z_p^2));

freqz(ans.Numerator{1}(:), ans.Denominator{1}(:));

%%

figure;
stem(1:length(hn),hn);
title('h[n]');
ylabel('Magnitude');
xlabel('n');
