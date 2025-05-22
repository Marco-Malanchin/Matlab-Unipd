clc
clear all
close all

N = 10;
T = 5;
Tp = N*T;

s1 = rand(1, N) + 1j*rand(1, N); %rand genera un vettore aleatorio con distribuzione standard

n_periodi_time = 3;

nT = [0:T:(n_periodi_time-1)*T];

figure;
subplot(2,1,1)
stem(nT,real(repmat(s1, 1, n_periodi_time)), "filled", "blu");
xlabel('nT');
title('parte reale')
subplot(2,1,2)
stem(nT,imag(repmat(s1, 1, n_periodi_time)), "filled", "red");
xlabel('nT');
title('parte immaginaria')

F = 1/Tp;
Fp = 1/T;

S1 = (1/N)*fft(s1); %trasformata di s1

n_periodi_freq = 2;

kF = [0:F:(n_periodi_freq*N-1)*F];

% figure;
% subplot(2,1,1)
% stem(kF, abs(repmat(S1, 1, n_periodi_freq)), 'black', 'Marker', '^');
% xlabel('kF [Hz]');
% title('Modulo Trasformata');
% 
% subplot(2,1,1)
% stem(kF, angle(repmat(S1, 1, n_periodi_freq)), 'black', 'Marker', '^');
% xlabel('kF [Hz]');
% title('Fase Trasformata');

s1_inv = N*ifft(S1);

s1_compare = [s1; s1_inv];

T = 10;
N = 20
Tp = N+T;

t = [0:T:(N-1)*T];

s2 = cos(2*pi*(1/tp)*t);

figure
stem(nT, repmat(s2, 1, n_periodi_time),'filled');
xlabel('nT [s]');

F = 1/Tp;
Fp = 1/t;

S2 =  (1/N)*fft(s2);

kF = [0:F:(n_periodi_freq)*F];

figure;
subplot(2,1,1);
stem(kF, abs(repmat(S2, 1, n_periodi_freq)),'black', 'Marker','^');
title('modulo trasformata');
subplot(2,1,1);
stem(kF, angle(repmat(S2, 1, n_periodi_freq)),'black', 'Marker','^');
title('fase trasformata');

kF_centrato = [-N*F:F(N-1)*F];
subplot(2,1,1);
stem(kF_centrato, abs(repmat(S2, 1, n_periodi_freq)),'black', 'Marker','^');
title('modulo trasformata');
subplot(2,1,1);
stem(kF_centrato, angle(repmat(S2, 1, n_periodi_freq)),'black', 'Marker','^');
title('fase trasformata');

%%
clc
clear all
close all

N = 10;
T = 5;
Tp = N*T;

s3 = rand(1,N) + 1j*rand(1, N);

S3_unitaria = (1/sqrt(N))*fft(s3);

%%
clc
clear all
close all

N = 1000;
T = 10;
Tp = N*T;

F = 1/Tp;

x = 2*rand(1, N) + 1j*(rand(1,N));

tic
X_dft = dft_manuale(x);
t_dft = toc;

tic
X_fft = (1/N)*fft(x);
t_fft = toc;

scarto = norm(X_dft - X_fft);

%%
N = 6;
a = [1; 2; 3; 4; 5; 6];