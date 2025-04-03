
function I = calcolatore_area_trapezio(f,t)

    if (length(f) == length(t))
        error("i vettori hanno lunghezze diverse!");
    end
        I = 0;
    for a = 1:length(f)-1
        I = I + (f(a+1) + f(a)) * (t(a+1) - t(a))/2;
    end

end
