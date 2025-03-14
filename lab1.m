a = 1;
b = 2;

c1 = a + b; %somma
c2 = a - b;
c3 = a * b;
c4 = a / b;
c5 = a^b;
c6 = sqrt(b);

disp(c1);
disp(c4);

stringa_prodotto = "il prodotto tra a e b è: %2.3f";
fprintf(stringa_prodotto, c3);

%%
clear all
clc

a = 1 + i;
b = 3 + 4i;
c = complex(10, 11);
d = exp(i*pi/4);
e = i^2;
f = sqrt(-1);

somma = a + b;
modulo_a = abs(a);
fase_a = angle(a);
coniugato_a = conj(a);
reale_d = real(d);
immaginario_d = imag(d);

%%
clear all 
clc

a = 57;
b = 30;

if (a == b)
    fprintf("a e b sono uguali");
else
    fprintf("a e b sono diversi");
end

if(a > b)
    fprintf(" a maggiore di b");
elseif (a < b)
    fprintf("a minore di b");
else
    fprintf("a uguale a b");
end

if (a > 10 & b > 10)
    fprintf("\n sia a che b sono maggiori di 10");
end

%%
clc
n = 5;
s = 0;
i = 1;
step = 1;

for z=1:step:n
    s = s+10;
    fprintf("\n iterazione n° %d, s = %2.1f \n", i, s);
    i = i+1;
end

%%
n = 5;
s = 0;
i = 1;
z = 1;
step = 1;

while (z <= n)
    s = s + 10;
    fprintf("Iterazione n° %d, s = %1.2f \n", i, s);
    i = i + 1;
    z = z + step;
end

%%
clc
clear all

v = [10, 20, 30, 40];
disp(v(1)); %il vettore parte da 1 non da 0
disp(v(end));
disp(v(2:end));

m = [10, 20, 30; 40, 50, 60; 70, 80, 90];
disp(m);
disp(m(1,3));

disp(m(1,:));
disp(m(:,3));

%%
clc
clear all

f = 5; %Hz
f2 = 10;
t_inizio = 0; %s
t_fine = 1; %s

fs = 22222220;
ts = 1 / fs;
t = [t_inizio:ts:t_fine]; %il vettore parte da tinizio e arrivo a tfine incrementando di ts ogni volta
y = sin(2*pi*f*t);
y2 = cos(2*pi*f2*t);

figure;
plot(t, y, '-');
xlabel("tempo [s]");
ylabel("ampiezza");
title("sinusoide");
grid on
hold on %sovrappone i grafici
plot(t, y2, '-');


