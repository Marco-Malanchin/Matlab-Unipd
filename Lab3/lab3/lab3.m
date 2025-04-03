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
f = k*exp(-alpha*t);

figure;
plot(t, f, 'Color', 'green');
xlabel('Tempo [s]');
ylabel('Ampiezza');
title('f(t)');
xlim([t_inizio-0.5 t_fine+0.5]);
ylim([min(f)-0.5 max(f)+0.5]);
grid on

A_trapezi = calcolatore_area_trapezio(f,t);
A = (-k/alpha)*(exp(-alpha*t_fine)) - exp(-alpha*t_inizio);

fprintf("La differenza (analitica-trapezi)è: %3.12f", A - A_trapezi);

A_sum = dt*sum(f);

fprintf("\nLAdifferenza (analitica - sum) è: %3.12f", A - A_sum);

%%
clc
clear all
close all

dt = 0.001;
t = [0:dt:3];

s = sin(2 * pi* t).*exp(-7 * (t-0.5).^2);

figure
plot(t, s, 'r--');
grid on

energia_s = dt*sum(abs(s).^2); %V^2*s

energia_s_norm = dt*(norm(s))^2;

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

s = 5*sin(2*pi*f1*t + pi/2) + 3*cos(2*pi*f2*t) + 2*sin(2*pi*f3*t);

figure;
plot(t, s);

periodo = 1/f1;

Ntp = floor(periodo/dt);

potenza_periodo = sum(abs(s(1:Ntp)).^2) / Ntp; 

%%
clc
clear all
close all

t_start = -5;
t_end = 5;
dt = 0.01;

time = [t_start:dt:t_end];

%primo rect
t11 = -0.8;
t21 = 0.8;
rect1 = zeros(1, length(time));
i_t11 = (t11 - t_start)/dt + 1;
i_t21 = (t21 - t_start)/dt +1;
rect1(i_t11:i_t21) = 1;


figure
plot(time, rect1);

t12 = -3;
t22 = 3;
rect2 = zeros(1, length(time));
i_t12 = (t12 - t_start)/dt + 1;
i_t22 = (t22 - t_start)/dt +1;
rect2(i_t11:i_t21) = 1;

figure;
plot(time, rect2);

[y, t_y] = calcolatore_convoluzione(rect1, time, rect2, time);

figure;
plot(t_y,y);

y_conv = dt*conv(rect1, rect2);

figure;
plot(t_y, y_conv);
