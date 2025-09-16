%% DEFINIZIONE SEGNALE
clc 
close all
clear all

N = 100;
T1 = 10;
dt = 2*T1/N;
t = [0:N-1]*dt-T1;

A = 2;
T2 = 5;

s = A*(abs(t) <= T2); %rettangolo discreto centrato in 0

%% RAPPRESENTAZIONE TEMPO

figure;
subplot(2,2,1)
plot(t,s)
xlabel('Tempo [s]');
ylabel('Ampiezza [V]');
ylim([-0.5 A+0.5]);
grid on

%% FORMA ANLITICA IN FREQUENZA
fa = [-1:0.001:1]; %Crea il vettore di frequenze fa, da –1 Hz a +1 Hz, con passo 0.001 Hz.
Sa = A*2*T2*sinc(fa*(2*T2)); %Questa è la trasformata di Fourier del rettangolo 

subplot(2,2,4)
plot(fa, real(Sa), 'color', 'blue', 'DisplayName', 'Parte Reale');
hold on
plot(fa, imag(Sa),'Color','red', 'DisplayName','Parte Immaginaria');
grid on
xlabel('Frequenza [Hz]')


%% USO fft SENZA FARE NULLA
S = fft(s) %trasformata del rettangolo senza aggiustare niente
subplot(2,2,2)
plot(real(S),'Color','blue','DisplayName','Parte reale');
hold on 
plot(imag(S),'Color','red','DisplayName','Parte immaginaria');
grid on 
xlabel('Frequenza [Hz]')
ylabel('[Vs]');
legend

%SOLUZIONE ESATTA
clc
clear all
close all

N = 1000; 
T1 = 70;
dt = 2*T1/N; 
t = [0:N-1]*dt-T1;

A = 2;
T2 = 5;
s = A*(abs(t) <= T2);

%Segnale nel tempo
figure;
subplot(2,2,1)
plot(t, s)
xlabel('Tempo [s]');
ylabel('Ampiezza [V]');
ylim([-0.5 A+0.5]);
grid on 

%Parte analitica
fa = [-3:0.001:3];

Sa = A*2*T2*sinc(fa*(2*T2)); %formula analitica del sinc

subplot(2,2,4)
plot(fa, real(Sa),'Color','blue','DisplayName','Parte reale');
hold on 
plot(fa, imag(Sa),'Color','red','DisplayName','Parte immaginaria');
grid on 
xlabel('Frequenza [Hz]')
xlim([-3 3])
ylabel('[Vs]');
legend
%

s = [s(N/2+1:N), s(1:N/2)]; %Stai riordinando il segnale nel tempo: prendi la seconda metà del vettore s e la metti davanti.

È lo stesso effetto di fftshift(s).

Risultato: il tuo rettangolo (che prima era definito da t=-T1 a t=+T1) viene riallineato in modo che l’indice 0 del vettore corrisponda al centro (t=0
t=0).

Così la rappresentazione nel tempo è centrata.

df = 1/(N*dt); %risoluzione in frequenza
f = [-N/2:N/2-1]*df; % crea il vettore delle frequenze (in Hz) ordinato da −𝐹𝑠/2 ad 𝐹𝑠/2−𝑑𝑓
👉 Questo serve per avere l’asse delle frequenze corretto quando poi disegni lo spettro.

S = dt*fft(s); %FFT del segnale rettangolare
S = [S(N/2+1:N), S(1:N/2)]; %centro lo spettro in 0

subplot(2,2,2)
plot(f, real(S),'Color','blue','DisplayName','Parte reale');
hold on 
plot(f, imag(S),'Color','red','DisplayName','Parte immaginaria');
grid on 
xlabel('Frequenza [Hz]')
ylabel('[Vs]');
xlim([-3 3])
legend


clc
clear all
close all

N = 120; 
T1 = 70;
dt = 2*T1/N; 
t = [0:N-1]*dt-T1;

A = 2;
T2 = 5;
s = A*(abs(t) <= T2);

%Segnale nel tempo
figure;
subplot(2,2,1)
plot(t, s)
xlabel('Tempo [s]');
ylabel('Ampiezza [V]');
ylim([-0.5 A+0.5]);
grid on 

%Parte analitica
fa = [-3:0.001:3];

Sa = A*2*T2*sinc(fa*(2*T2));

subplot(2,2,4)
plot(fa, real(Sa),'Color','blue','DisplayName','Parte reale');
hold on 
plot(fa, imag(Sa),'Color','red','DisplayName','Parte immaginaria');
grid on 
xlabel('Frequenza [Hz]')
xlim([-3 3])
ylabel('[Vs]');
legend
%

s = [s(N/2+1:N), s(1:N/2)];

df = 1/(N*dt);
f = [-N/2:N/2-1]*df;

S = dt*fft(s);
S = [S(N/2+1:N), S(1:N/2)];

subplot(2,2,2)
plot(f, real(S),'Color','blue','DisplayName','Parte reale');
hold on 
plot(f, imag(S),'Color','red','DisplayName','Parte immaginaria');
grid on 
xlabel('Frequenza [Hz]')
ylabel('[Vs]');
xlim([-3 3])
legend
