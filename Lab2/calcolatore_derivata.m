function fd = calcolatore_derivata(f,t)
    fd = zeros(1, length(f));
    for n = 1:length(t)-1
        fd(n) = (f(n+1)-f(n))/(t(n+1)-t(n));
    end

    fd(end) = fd(end) - 1;%l'ultimo elemento del vettore sarà uguale al penultimo elemento calcolato, perchè nel for calcolo usando il successivo, ma l'ultimo elemento non avrà un sucessivo
end

