% %SECONDO METODO
% fd = diff(f) ./ diff(t);
%  fd = [fd, fd(end)]; 

%area di un segnale
clc
clear all
close all

t_inizio = 0;
t_fine = 10;
dt = 0.05;
t = [t_inizio:dt:t_fine]; %s

alpha = 0.1;
k = 3;
f = k*exp(-alpha*t); %segnale esponenziale smorzato

figure;
plot(t, f, 'Color', 'green');
xlabel('Tempo [s]');
ylabel('Ampiezza');
title('f(t)');
xlim([t_inizio-0.5 t_fine+0.5]);
ylim([min(f)-0.5 max(f)+0.5]);
grid on

A_trapezi = calcolatore_area_trapezio(f,t);
A = (-k/alpha)*(exp(-alpha*t_fine) - exp(-alpha*t_inizio));

fprintf("La differenza (analitica-trapezi)è: %3.12f", A - A_trapezi); %Stampa l’errore tra integrale analitico e quello numerico a trapezi. Se la frequenza di campionamento fs è alta, la differenza sarà molto piccola.
% --- Metodo delle somme rettangolari ---
% Approssimazione dell'integrale come somma di rettangoli
A_sum = dt*sum(f);

fprintf("\nLAdifferenza (analitica - sum) è: %3.12f", A - A_sum);

%%
clc
clear all
close all

dt = 0.001;
t = [0:dt:3];

s = sin(2 * pi* t).*exp(-7 * (t-0.5).^2); %Segnale definito come un sinusoide modulata da una gaussiana centrata in t=0.5.

figure
plot(t, s, 'r--');
grid on

energia_s = dt*sum(abs(s).^2); %V^2*s %Calcola l’energia del segnale con definizione discreta, abs(s).^2 → valore al quadrato del segnaledt*sum(...) → moltiplica per l’intervallo di campionamento, cioè approssima l’integrale.

energia_s_norm = dt*(norm(s))^2; %Qui calcoli la stessa energia usando la norma 2 di MATLAB.

%
clc
clear all
close all

dt = 0.001;
t_inizio = 0;
t_fine = 3;
t = [t_inizio:dt:t_fine];

f1 = 5;
f2 = f1*2;
f3 = f1*3;

s = 5*sin(2*pi*f1*t + pi/2) + 3*cos(2*pi*f2*t) + 2*sin(2*pi*f3*t);   % prima armonica, ampiezza 5, sfasata di 90°. seconda armonica, ampiezza 3. terza armonica, ampiezza 2

figure;
plot(t, s);

periodo = 1/f1; % periodo fondamentale = 0.2 s

Ntp = floor(periodo/dt);  % numero di campioni in un periodo

potenza_periodo = sum(abs(s(1:Ntp)).^2) / Ntp; %Calcolo della potenza media del segnale usando un solo periodo:

%%
clc
clear all
close all

t_start = -5;
t_end = 5;
dt = 0.01;

time = [t_start:dt:t_end]; % vettore tempo (1001 campioni)

% --- Primo rettangolo: da -0.8 a 0.8 ---
t11   = -0.8;
t21   =  0.8;
rect1 = zeros(1, length(time));
i_t11 = round((t11 - t_start)/dt) + 1;   % indice corrispondente a t11
i_t21 = round((t21 - t_start)/dt) + 1;   % indice corrispondente a t21
rect1(i_t11:i_t21) = 1;


figure
plot(time, rect1);

% --- Secondo rettangolo: da -3 a 3 ---
t12   = -3;
t22   =  3;
rect2 = zeros(1, length(time));
i_t12 = round((t12 - t_start)/dt) + 1;
i_t22 = round((t22 - t_start)/dt) + 1;
rect2(i_t12:i_t22) = 1; 

figure;
plot(time, rect2);

[y, t_y] = calcolatore_convoluzione(rect1, time, rect2, time);

figure;
plot(t_y,y);

% --- Convoluzione "manuale" con conv + asse tempi corretto ---
y_conv  = dt * conv(rect1, rect2);%conv Fa la convoluzione discreta tra due vettori. in analisi dei segnali, la convoluzione è un integrale Quando la implementi con somme discrete (conv), bisogna moltiplicare per dt per approssimare correttamente l’integrale:
t_conv  = (t_start + t_start) : dt : (t_end + t_end); % da -10 a 10, passo dt

figure
plot(t_conv, y_conv, 'LineWidth', 1.5);
grid on
title('Convoluzione (dt * conv)')
