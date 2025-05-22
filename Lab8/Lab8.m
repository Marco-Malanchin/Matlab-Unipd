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

s = A*(abs(t) <= T2);

%% RAPPRESENTAZIONE TEMPO

figure;
subplot(2,2,1)
plot(t,s)
xlabel('Tempo [s]');
ylabel('Ampiezza [V]');
ylim([-0.5 A+0.5]);
grid on

%% FORMA ANLITICA IN FREQUENZA
fa = [-1:0.001:1];
Sa = A*2*T2*sinc(fa*(2*T2));

subplot(2,2,4)
plot(fa, real(Sa), 'color', 'blue', 'DisplayName', 'Parte Reale');
hold on
plot(fa, imag(Sa),'Color','red', 'DisplayName','Parte Immaginaria');
grid on
xlabel('Frequenza [Hz]')


%% USO fft SENZA FARE NULLA
S = fft(s)