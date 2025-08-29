%%
clc;
clear;
close all; %chiude tutte le figure

v = [10,20,30,40];
disp(v);
disp(v(3));
disp(v(end));
disp(v(end-1));

%%
clc
clear
close all
m = [10, 20, 30; 40, 50 ,60; 70, 80, 90; 100, 110, 120]; %matrice 4x3
dismp(m);
disp(m(3,2)); % stampa l’elemento della terza riga, seconda colonna di m
              % cioè 80
disp(m(3,:));
disp(m(:,1));
disp(size(m));
disp(size(m,3));

%%
clc
clear
close all
n = zeros(2,3); %matrice 2x3 fatta di 0
disp(n);
n1 = ones(2,3);
disp(n1);
n2 = eye(5); %crea una matrice identità 5×5 e la assegna alla variabile n2.
disp(n2);

%%
clc
clear
close all

t_inizio = 0; %secondi
t_fine = 1; %secondi
fs = 1000; %frequenza di campionamento
ts = 1/fs;
t = [t_inizio:ts:t_fine]; %crea un vettore che parte da 0, cresce a passi di 0.001 secondi e arriva fino a 1 secondo. In totale hai (t_fine - t_inizio)/ts + 1 = (1 - 0)/0.001 + 1 = 1001 campioni.
F = 5; % frequenza della sinusoide
y = sin(2*pi*F*t);
y2 = cos(2*pi*F*t);
figure;
plot(t,y,'g', LineWidth= 3); %'g' e' il colore che uso per plottare il grafico, lineWidrth aumenta lo spessore della linea del grafico
xlabel('tempo [s]');
ylabel('ampiezza');
ylim([-2 2]); %aumento la larghezza della visione dell asse y
xlim([-1 2]); %aumento la larghezza della visione dell asse x
grid on
hold on
plot(t,y2,'b', LineWidth= 3);%plotto il grafico del coseno sullo stesso grafico del seno

%%
clc
clear
close all

delta_t = 0.5e-6; %tempo che passa tra l'emissione e la ricezione di un segnale es radar.
distanza = calcolatore_distanza(delta_t); % funzione che abbiamo creato per calcolare la distenza in base al delta temporale
fprintf("La distanza è %2.4f metri\n", distanza);

%%
clc
clear
close all

t_start = 0;
t_end = 6;
t_step = 0.001;

t = [t_start:t_step:t_end];
F = 1;
f = sin(2*pi*F*t);
fd = calcolatore_derivata(f,t);

figure;
plot(t,f);
hold on
plot(t, fd);
