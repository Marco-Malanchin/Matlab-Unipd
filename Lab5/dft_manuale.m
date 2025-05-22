function X_dft = dft_manuale(x)
    N = length(x);

    X_dft = zeros(1, N);
    for k = 0:N-1
        for n = 0:N-1
            X_dft(k+1) = X_dft(k+1) + (1/N) * x(n+1) * exp(-1j*2*pi*k*n/N);

        end
    end
end