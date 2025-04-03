

function [y, t_y] = calcolatore_convoluzione(x, t_x, h, t_h)
    
    dt_x = t_x(2)-t_x(1);
    dt_h = t_h(2)-t_h(1);

    Lx = length(x);
    Lh = length(h);
    Ly = Lx + Lh -1;

    t_min = t_x(1) + t_h(1);

    t_y = t_min + (0:Ly-1)*dt;

    y = zeros(1, Ly);

    for n = 1:Ly
        for k = 1:Lh
            if (n-k+1 > 0) && (n-k+1 <= Lx)
                y(n) = y(n) + (dt*x(n-k+1)*h(k));
            end

        end
    end
end