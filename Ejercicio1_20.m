format long

% Funciones de calculo de la TIR.
% Datos iniciales: intereses y tiempo maximo.

y=[0.03,0.03,0.03,0.03,0.03,0.0325,0.035,0.0375,0.04,0.0425];
tm=10;

% Precio de emision de los bonos, caso compuesto y continuo,
% respectivamente
B_1=B_comp(y,tm);
B_2=B_cont(y,tm);


% y de la TIR en el caso compuesto
fprintf("Valor de y para el caso compuesto")
Y1=dicot_TIR(0.03,0.05,@TIR_comp,B_1,tm)
% Comprobacion para el caso continuo
p1=TIR_comp(Y1,B_1,tm);

der_B1=d_TIR_comp(Y1,tm);
der2_B1=d2_TIR_comp(Y1,tm);

% Duracion para el caso compuesto
fprintf("Duracion para el caso compuesto")
D1=-(1+Y1)*der_B1/B_1

% Convexidad para el caso compuesto
fprintf("Convexidad para el caso compuesto")
C1=der2_B1/B_1



% y de la TIR en el caso continuo
fprintf("Valor de y para el caso continuo")
Y2=dicot_TIR(0.03,0.05,@TIR_cont,B_2,tm)
% Comprobacion para el caso continuo
p2=TIR_cont(Y2,B_2,tm);

der_B2=d_TIR_cont(Y2,tm);
der2_B2=d2_TIR_cont(Y2,tm);

% Duracion para el caso compuesto
fprintf("Duracion para el caso continuo")
D1=-(1+Y2)*der_B2/B_2

% Convexidad para el caso compuesto
fprintf("Convexidad para el caso continuo")
C1=der2_B2/B_2


function c=dicot_TIR(a,b,TIR,B_c,tm)
    c=(a+b)/2;
    fa=TIR(a,B_c,tm);
    fb=TIR(b,B_c,tm);
    fc=TIR(c,B_c,tm);
    if fa*fb<0
        %fprintf("Hay un cambio de signo en el intervalo [a,b]. Aplicando dicotomia a la funcion.")
        for i=1:100         % Suponemos por defecto 100 iteraciones.
            if fa*fc<0
                b=c;
                c=(a+b)/2;
                fb=TIR(b,B_c,tm);
                fc=TIR(c,B_c,tm);
            elseif fc*fb<0
                a=c;
                c=(a+b)/2;
                fa=TIR(a,B_c,tm);
                fc=TIR(c,B_c,tm);
            end
        end
    else
        fprintf("No hay cambio de signo en la funcion. Devolvemos el valor 0 por defecto.")
        c=0;
    end    
end


% FUNCIONES PARTE COMPUESTA

% Funcion precio de emision de bono compuesta
function f=B_comp(y,tm)
    f=0;
    for i=1:tm-1
        f=f+6/(1+y(i))^i;
    end
    f=f+1006/(1+y(tm))^tm;
end

% Funcion que resta la TIR con el valor de emision del bono en el caso compuesto, para poder
% aplicar el metodo de dicotomia.
function f=TIR_comp(y,B,tm)
    f=B;
    for i=1:tm-1
        f=f-6/(1+y)^i;
    end
    f=f-1006/(1+y)^tm;
end

%derivada primera de la TIR compuesta
function f=d_TIR_comp(y,tm)
    f=0;
    for i=1:tm-1
        f=-f+6*i/(1+y)^(i+1);
    end
    f=f-1006*tm/(1+y)^(tm+1);
end

%derivada segunda de la TIR compuesta

function f=d2_TIR_comp(y,tm)
    f=0;
    for i=1:tm-1
        f=f+6*i*(i+1)/(1+y)^(i+2);
    end
    f=f+1006*tm*(tm+1)/(1+y)^(tm+2);
end

% FUNCIONES PARTE CONTINUA

% Funcion precio de emision de bono continuna
function f=B_cont(y,tm)
    f=0;
    for i=1:tm-1
        f=f+6*exp(-y(i)*i);
    end
    f=f+1006*exp(-y(tm)*tm);
end

% Funcion que resta la TIR con el valor de emision del bono en el caso continuo, para poder
% aplicar el metodo de dicotomia.
function f=TIR_cont(y,B,tm)
    f=B;
    for i=1:tm-1
        f=f-6*exp(-y*i);
    end
    f=f-1006*exp(-y*tm);
end

% Funcion derivada primera de la TIR continua
function f=d_TIR_cont(y,tm)
    f=0;
    for i=1:tm-1
        f=f-6*exp(-y*i)*i;
    end
    f=f-1006*exp(-y*tm)*tm;
end

% Funcion derivada segunda de la TIR continua
function f=d2_TIR_cont(y,tm)
    f=0;
    for i=1:tm-1
        f=f+6*exp(-y*i)*i^2;
    end
    f=f+1006*exp(-y*tm)*tm^2;
end