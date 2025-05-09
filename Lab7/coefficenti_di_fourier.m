clc
clear all
close all
N = 30;
Tp = 1;

dt = Tp/N;

t = [0:N-1]*dt;

fcos = 1;
s = cos(2*pi*fcos*t);

figure;
plot(t, s);
xlabel('t [s]')
ylabel('ampiezza')
grid on

S = dt/Tp*fft(s);
S = [S(N/2 + 1:N), S(1:N/2)]; %creo un vettore con la trasformata appena creata, mettendo la seconda metà come primo elemento e la prima metà come secondo, in questo modo la trasformata sarà graficamente giusta
F = 1/Tp;
f = [-N/2:N/2-1]*F;
figure;
stem(f, abs(S), 'filled');
xlabel('Frequenza [Hz]')
grid on
ylim([-0.5 1])

%%
clc
clear all
close all
N = 500;
Tp = 1;

dt = Tp/N;

t = [0:N-1]*dt;

s = zeros(1, N);
s(1:N/2) = 1;
s(N/2+1:N) = 0;

figure
plot(t, s)
ylim([-0.5 1.5])
xlabel('t[s]')
ylabel('ampiezza')
grid on

S = dt/Tp*fft(s);
S = [S(N/2+1:N), S(1:N/2)];

F = 1/Tp;
f = [-N/2:N/2-1]*F;

figure
stem(f, abs(S), 'filled');
xlabel('Frequenza [Hz]')
grid on


