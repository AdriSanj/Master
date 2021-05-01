%%%% Programa para calcular la TAE.

% Lista de tipos de acumulacion. Elegir en la linea 8 la posicion del tipo
% de acumulacion en la lista.
acumulacion=["anual","semestral","trimestral","mensual","semanal","diaria","horaria"];

fprintf("Acumulacion compuesta")
TAE_comp(2.5,acumulacion(1))        %% Llamada a la funcion que calcula la TAE compuesto, con un interes de 2.5%

frpintf("Acumulacion continua")     %% Llamada a la funcion que calcula la TAE continua, con un interes de 2.5%
TAE_cont(2.5)

function yp=TAE_comp(y,acumulacion)
    y=y/100;
    if acumulacion == "anual"
        m=1;
    elseif acumulacion == "semestral"
        m=2;
    elseif acumulacion == "trimestral"
        m=4;
    elseif acumulacion == "mensual"
        m=12;
    elseif acumulacion == "semanal"
        m=52;
    elseif acumulacion == "diaria"
        m=365;
    elseif acumulacion == "horaria"
        m=8760;
    end
    yp=(1+y/m)^m-1;
end

function yp=TAE_cont(y)
    y=y/100;
    yp=exp(y)-1;
end