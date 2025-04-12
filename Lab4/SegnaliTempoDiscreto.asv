%%
clc
clear all
close all
T = 2; %s
nT = [-10:T:50];%vettore campionamenti
k = 5; %traslazione
nT_traslato = nT-k;

alpha = 0.9;%coefficente
f = 1/12; %Hz
x = (nT >= 0) .* ((alpha.^nT)) .* exp(1j * 2 * pi * f *nT); %il punto "." prima delle opoerazioni, applica l'operazione ad ogni elemento
x_traslato= (nT_traslato >= 0) .* ((alpha.^nT_traslato)) .* exp(1j * 2 * pi * f *nT_traslato)
x_reale = real(x);
x_immaginario = imag(x);
x_modulo = abs(x);
x_fase = angle(x);

figure;
subplot(2,2,1);
stem(nT, x_reale,"filled", 'DisplayName', 'Parte reale'); %comando per rappresnetare i segnali a tempo discreto
grid on
hold on
stem(nT_traslato, real(x_traslato),"filled", 'DisplayName', 'Parte reale');
xlabel('nT [secondi]');
ylabel('Ampiezza');
legend;


subplot(2,2,2);
stem(nT, x_immaginario,"filled", 'DisplayName', 'Parte immaginaria'); %comando per rappresnetare i segnali a tempo discreto
grid on
hold on
stem(nT_traslato, imag(x_traslato),"filled", 'DisplayName', 'Parte immaginaria');
xlabel('nT [secondi]');
ylabel('Ampiezza');
legend;


subplot(2,2,3);
stem(nT, x_modulo,"filled", 'DisplayName', 'Modulo'); %comando per rappresnetare i segnali a tempo discreto
grid on
hold on
stem(nT_traslato, abs(x_traslato),"filled", 'DisplayName', 'Modulo');
xlabel('nT [secondi]');
ylabel('Ampiezza');
legend;


subplot(2,2,4);
stem(nT, x_fase,"filled", 'DisplayName', 'Fase'); %comando per rappresnetare i segnali a tempo discreto
grid on
hold on
stem(nT_traslato, angle(x_traslato),"filled", 'DisplayName', 'Fase');
xlabel('nT [secondi]');
legend;

energia_segnale = T*sum(abs(x).^2);

%%
clc
clear all
close all

T = 1;
nT = [-20:T:20]; %s

f = 1/20; %Hz
x = cos(2*pi*f*nT);

figure;
stem(nT,x, "filled");
grid on;
xlabel('Tempo[s]');
ylabel('Ampiezza');

Tp = 1/f;
campioni_per_periodo = floor(Tp/T);
E_priodo = T*sum(abs(x(1:campioni_per_periodo)).^2);

%%
clc
clear all
close all

T = 1;
Np = 10;

t = T*(-floor(Np/2) : floor((Np-1))/2);
D = 0.5;
L_rect = round(D*Np);
x = zeros(1, Np);
start_idx = ceil((Np-L_rect)/2)+1;
x(start_idx : start_idx + L_rect - 1) = 1;

h = x;

figure;
stem(t, h, 'filled');
xlabel('nT [s]');
grid on;
ylim([0 1.2]);

y = T*conv(x, h);
Ly = length(y);

t_y = (2*t(1) + T*(0:Ly-1));

figure;
stem(t_y, y, 'filled');
xlabel('nT [s]');
grid on;
title('convoluzione');

N_periodi = 7;
x_p = repmat(x, 1, N_periodi);
h_p = repmat(h, 1, N_periodi);

t_p = [];
for k = -floor(N_periodi/2) : floor((N_periodi-1)/2)
    t_p = [t_p, t+ k*Np];
end

figure;
stem(t_p, x_p, 'filled');
xlabel('nT [s]');
ylabel('Ampiezza');
ylim([0 1.2]);
grid on;

figure;
stem(t_p, h_p, 'filled');
xlabel('nT [s]');
ylabel('Ampiezza');
ylim([0 1.2]);
grid on;

y_p = T*conv(x_p, h_p);
Ly_p = length(y_p);
t_y_p = t_p(1) + t_p(1) + T*(0:Ly_p-1);

figure;
stem(t_y_p, y_p, 'filled');
xlabel('nT [s]');
ylabel('Ampiezza');
grid on;
title('Convoluzione dei Segnali Periodici')

h_rep2 = repmat(h,1,2);
y_rep2 = T*conv(x, h_rep2);

t_h_rep2 = [];
for k = -floor(2/2) : floor((2-1)/2)
    t_h_rep2 = [t_h_rep2, t + k*Np];
end

t_y_rep2 = t(1) + t_h_rep2(1) + T*[0: length(y_rep2)-1];

figure;
stem(t_y_rep2, y_rep2, 'filled');


