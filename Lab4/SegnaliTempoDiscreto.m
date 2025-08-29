%%
clc
clear all
close all
T = 2; % [secondi] periodo di campionamento (passo tra i campioni)
nT = [-10:T:50]; %vettore campionamenti
k = 5; %traslazione temporale (ritardo positivo)
nT_traslato = nT-k; %griglia temporale traslata

alpha = 0.9; %coefficiente di decadimento (0<alpha<1 => decresce per t>=0)
f = 1/12; %Hz
% Segnale: x(t) = u(t) * alpha^t * e^{j 2π f t}
% u(t) = (t>=0) impone causalità: x(t<0)=0
x = (nT >= 0) .* ((alpha.^nT)) .* exp(1j * 2 * pi * f *nT); %il punto "." prima delle opoerazioni, applica l'operazione ad ogni elemento
% Segnale traslato: x(t-k) = u(t-k) * alpha^{t-k} * e^{j 2π f (t-k)}
% ATTENZIONE: qui stai valutando la formula su una griglia TEMPO TRASLATA,
% non stai "shiftando di m campioni" il vettore. È una traslazione continua di k secondi.
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

energia_segnale = T*sum(abs(x).^2); %formula energia segnale tempo discreto

%%
clc
clear all
close all

T = 1; % passo di campionamento (1 s)
nT = [-20:T:20]; %istanti di campionamento [-20, -19, ..., 20]

f = 1/20; %Hz
x = cos(2*pi*f*nT); %segnale discreto campionato

figure;
stem(nT,x, "filled");
grid on;
xlabel('Tempo[s]');
ylabel('Ampiezza');

Tp = 1/f; % periodo fondamentale (20 s)
campioni_per_periodo = floor(Tp/T); % numero di campioni in un periodo 20
E_priodo = T*sum(abs(x(1:campioni_per_periodo)).^2); %prendi i primi 20 campioni del vettore x.
→ Siccome nT va da -20 a 20, i “primi 20” non corrispondono esattamente a un periodo intero, ma a un intervallo qualsiasi.

Però per un coseno periodico non cambia, perché ogni periodo ha la stessa energia.

%%
clc
clear all
close all

T = 1;
Np = 10;

t = T*(-floor(Np/2) : floor((Np-1))/2); %con Np=10 -> t = -5:4 (10 campioni)
D = 0.5; % duty cycle = 50% -> il rect dura metà del periodo
L_rect = round(D*Np); % lunghezza (in campioni) del rect = 0.5*10 = 5, Vuol dire che il tuo periodo è fatto da 10 campioni (Np=10), e il rettangolo deve durare 5 campioni (cioè metà periodo).
x = zeros(1, Np); %Vettore di 10 campioni inizialmente tutti a zero.
start_idx = ceil((Np-L_rect)/2)+1;
%Np-L_rect = 10-5 = 5 → “quanto spazio libero resta da riempire con zeri”.

%(Np-L_rect)/2 = 2.5 → metà dello spazio libero, cioè per centrare il rect.

%ceil(...)+1 → ti dà l’indice iniziale dove cominciare a piazzare gli 1.
%In questo caso: 4 (quarto campione del vettore).
x(start_idx : start_idx + L_rect - 1) = 1;
%Parti dall’indice 4 e scrivi L_rect=5 valori a 1.

%Quindi metti a 1 i campioni 4,5,6,7,8.

%👉 Il tuo x finale sarà:

[0 0 0 1 1 1 1 1 0 0]


%cioè un rettangolo di 5 campioni centrato in un periodo da 10 campioni.

h = x;

figure;
stem(t, h, 'filled');
xlabel('nT [s]');
grid on;
ylim([0 1.2]);

y = T*conv(x, h); %scala per T per interpretazione “integrale” conv discreta
Ly = length(y);  

t_y = (2*t(1) + T*(0:Ly-1)); % tempo parte da t_start_x + t_start_h = -10

figure;
stem(t_y, y, 'filled');
xlabel('nT [s]');
grid on;
title('convoluzione');

N_periodi = 7; % Quante copie consecutive dei segnali vuoi
x_p = repmat(x, 1, N_periodi); % x ripetuto N_periodi volte (concatenazione)
h_p = repmat(h, 1, N_periodi); % h ripetuto N_periodi volte
%repmat(...,1,N) concatena i vettori N volte → ottieni due sequenze “lunghe”.
t_p = [];
for k = -floor(N_periodi/2) : floor((N_periodi-1)/2)
    t_p = [t_p, t+ k*Np];
end
%N_periodi → numero di copie (periodi) che vuoi replicare di 
%x e h.
%t → asse temporale di un solo periodo (lunghezza Np campioni).
%Ciclo su k: ad ogni iterazione, trasli t di k*Np campioni e lo concateni.
%Alla fine ottieni un vettore t_p lungo Np*N_periodi, che rappresenta l’asse dei tempi per la sequenza periodica.

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
% --- Convoluzione lineare delle sequenze replicate ---
y_p = T*conv(x_p, h_p);
Ly_p = length(y_p);
% L’asse dei tempi della conv lineare parte dalla somma dei primi istanti:
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


