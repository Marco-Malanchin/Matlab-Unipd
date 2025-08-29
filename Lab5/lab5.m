clc
clear all 
close all

% Parametri tempo/frequenza
N  = 10;                 % campioni per periodo
T  = 5;                  % passo temporale [s]
Tp = N*T;                % durata periodo [s]

% Segnale periodico base (complesso)
s1 = rand(1, N) + 1j*rand(1, N);

n_periodi_time = 3;

nT = [0:T:(n_periodi_time*N-1)*T]; %nT diventa un vettore di tempi che contiene tutti i campioni successivi dei periodi replicati.

figure;
subplot(2,1,1)
stem(nT, real(repmat(s1, 1, n_periodi_time)), 'filled','blu');
grid on
xlabel('nT [s]')
title('Parte reale')
subplot(2,1,2)
stem(nT, imag(repmat(s1, 1, n_periodi_time)), 'filled','red');
grid on
xlabel('nT [s]')
title('Parte immaginaria')

F = 1/Tp;
Fp = 1/T;

S1 = (1/N)*fft(s1); %comando per la trasformata, manca 1/N al comando, dobbiamo aggiungerlo noi
% Repliche dello spettro (solo per visualizzazione)
n_periodi_freq = 2;

kF = [0:F:(n_periodi_freq*N-1)*F]; %Il vettore kF rappresenta l’asse delle frequenze corrispondente ai punti di S1_rep.

figure;
subplot(2,1,1)
stem(kF, abs(repmat(S1, 1, n_periodi_freq)),'black','Marker','^');
grid on
xlabel('kF [Hz]')
title('Modulo')
subplot(2,1,2)
stem(kF, angle(repmat(S1, 1, n_periodi_freq)),'magenta','Marker','^');
grid on
xlabel('kF [Hz]')
title('Fase')

s1_inv = N*ifft(S1); %inverso della trasformata fft, dobbiamo aggiungere N

s1_compare = [s1; s1_inv]; %Crea una matrice 2×N:

%Prima riga = segnale originale s1

%Seconda riga = ricostruzione da FFT+IFFT


figure;
subplot(2,1,1)
stem(nT, real(repmat(s1_inv, 1, n_periodi_time)), 'filled','blu');
grid on
xlabel('nT [s]')
title('Parte reale inversa')
subplot(2,1,2)
stem(nT, imag(repmat(s1_inv, 1, n_periodi_time)), 'filled','red');
grid on
xlabel('nT [s]')
title('Parte immaginaria inversa')

T = 10;
N = 20;
Tp = T*N;

t = [0:T:(N-1)*T]; %se ho 20 campioni e il primo parte da 0, significa che l'ultimo sarà N-1

s2 = cos(2*pi*(1/Tp)*t);

nT = [0:T:(n_periodi_time*N-1)*T];

figure;
stem(nT, repmat(s2, 1, n_periodi_time), 'filled');
xlabel('nT [s]');
grid on
F = 1/Tp;
Fp = 1/T;

S2 = (1/N)*fft(s2);

kF = [0:F:(n_periodi_freq*N-1)*F];

figure;
subplot(2,1,1)
stem(kF, abs(repmat(S2, 1, n_periodi_freq)),'black','Marker','^');
grid on
xlabel('kF [Hz]')
title('Modulo')
subplot(2,1,2)
stem(kF, angle(repmat(S2, 1, n_periodi_freq)),'magenta','Marker','^');
grid on
xlabel('kF [Hz]')
title('Fase')

kF_centered = [-N*F:F:(N-1)*F]; %Facciamo lo spettro “centrato” perché la FFT standard di MATLAB mette le frequenze negative in fondo, mentre nel dominio dei segnali è molto più comodo visualizzare da –Fs/2 a +Fs/2 attorno a 0 Hz.

figure;
subplot(2,1,1)
stem(kF_centered, abs(repmat(S2, 1, n_periodi_freq)),'black','Marker','^');
grid on
xlabel('kF [Hz]')
title('Modulo')
subplot(2,1,2)
stem(kF_centered, angle(repmat(S2, 1, n_periodi_freq)),'magenta','Marker','^');
grid on
xlabel('kF [Hz]')
title('Fase')

%%
clc
clear all
close all

N = 10;
T = 5;
Tp = N*T;

s1 = rand(1, N) + 1j*rand(1, N);

S1_unit = (1/sqrt(N))*fft(s1);

fprintf("Norma della trasformata unitaria: %3.3f", norm(S1_unit));
fprintf("\nNorma del segnale s1: %3.3f", norm(s1));

clc
clear all
close all

N = 1000;  % Lunghezza del segnale
T = 10;
Tp = N*T;

F = 1/Tp;

x = 2*rand(1, N)-1 + 1j*(2*rand(1, N)-1);


%manuale
tic; %->Avvia cronometro
X_dft = dft_manuale(x); %la si può implementare a mano con due cicli for, fft matlab molto più veloce, manuale(N^2) fft(N logN)
t_dft = toc;  %->Tempo DFT manuale


% fft
tic;%->Avvia cronometro
X_fft = (1/N)*fft(x);
t_fft = toc;  %->Tempo fft 

fprintf('Tempo FFT         : %.4f secondi\n', t_fft);

errore = norm(X_dft - X_fft);  % Norm L2
fprintf('Errore tra dft manuale e fft: %.2e\n', errore);

%%
N = 6;
a = [1; 2; 3; 4; 5; 6];
