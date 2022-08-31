%% 5)

ts=6400;
alfa=-cos((0.25*pi+0.7*pi)/2)/cos((0.7*pi-0.25*pi)/2);
Z=tf([-alfa -1],[1 alfa],ts,'variable','z^-1');

A=0.001836*(1+Z)^4;
B=(1-1.5548*Z+0.6493*Z^2)*(1-1.4996*Z+0.8482*Z^2);

C=A/B

%figure;
%freqz(C.Numerator{1}(:),C.Denominator{1}(:));

figure;
[H1 W1]=freqz(C.Numerator{1}(:),C.Denominator{1}(:));
plot(W1/pi,abs(H1));
ylabel('Magnitude');
xlabel('Frequency');
title('Highpass Filter');

%% 6)
k=1.2748;
alfa2=0.3249;
Z2=tf([-0.1208 0.3641 -1],[1 -0.3641 0.1208],ts,'variable','z^-1');

A2=0.001836*(1+Z2)^4;
B2=(1-1.5548*Z2+0.6493*Z2^2)*(1-1.4996*Z2+0.8482*Z2^2);

C2=A2/B2

%figure;
%freqz(C2.Numerator{1}(:),C2.Denominator{1}(:));

figure;
[H1 W1]=freqz(C2.Numerator{1}(:),C2.Denominator{1}(:));
plot(W1/pi,abs(H1));
ylabel('Magnitude');
xlabel('Frequency');
title('Bandpass Filter');

%% 7) 

delta=10^-5;
alfa=42.5;

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

figure
stem(x);
title("Old Plot");
xlabel('n');
ylabel('Magnitude');



x_new= x;
fact = 4;
pad_length = length(x_new)*(fact - 1);
z_length = ceil((length(x_new)+1)/2);
z = fft(x_new);

z_p = [z(1:z_length) zeros(1, pad_length) z(z_length+1:end)];
if  (mod(length(x_new), 2)==0)
    z_p(z_length) = z_p(z_length)/2; 
    z_p(z_length+pad_length) = z_p(z_length);
end

x_p = real(ifft(z_p))*fact; 

figure
stem(x_p);
title("New Plot");
xlabel('n');
ylabel('Magnitude');
