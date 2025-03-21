function d = calcolatore_distanza(delta_t)
    c = 3e8; %assegno a c il valore della velocità della luce
    if(d < 0)
        error('un tempo non puo essere negativo')
    else
    d = (delta_t*c)/2;%prendiamo il tempo lo moltiplichiamo per la velocità della luce e lo dividiamo per 2 perche contiamo sia andata che ritorno
end

