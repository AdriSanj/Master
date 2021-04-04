% Instrucciones de uso:
% El programa toma inicialmente las constantes del problema planteado. Se
% pueden modificar al gusto. En cada una viene indicado que magnitud fisica
% representa

% Ecuacion de ondas

c0=1500;            % Velocidad de propagacion de la onda en el medio
rho=1000;            % Densidad del medio

%Condiciones iniciales
ur=0;       % Velocidad a la derecha
ul=0;       % Velocidad a la izquierda
pl=1e6;     % Presion a la izquierda
pr=1e3;     % Presion a la derecha
%P=max([pl,pr]);        % En caso de querer otro tipo de escalado,
%descomentar esta linea junto con las lineas en las cuales los comandos
%axis estan comentados.

%Variables caracteristicas en el centro
alf1R=0.5*(ur-pr/(c0*rho));
alf2L=0.5*(ul+pl/(c0*rho));

N=1000;                  % Numero de puntos del mallado
x=linspace(-5,5,N);     % Vector de mallado
t=0;                    % Instante de tiempo inicial
mu=1          
Dx=(x(N)-x(1))/N;       % Espacio entre puntos del mallado
Dt=Dx*mu/c0;  % Espacio entre tiempos


% Ejes modificados de las graficas. Descomentar si se quiere para tener
% otro escalado, junto con lo que ya se ha hablado anteriormente

% xmin=x(1);
% xmax=x(N);
% yminp=-2*P;
% ymaxp=2*P;

% Eleccion del apartado que queremos representar. Posibilidades:
% apartado="c"
% apartado="d"
% apartado= "er" se corresponde con el apartado e y problema de riemann
% apartado= "ec" apartado e y problema de cauchy del apartdo d
% Es importante la eleccion de estas y solo estas opciones para el uso
% correcto del programa. En caso contrario no ejecuta los calculos.

apartado="er";
% Para el apartado e se pueden comparar los resultados de los apartados
% anteriores con los de este, eligiendo la opcion "si". En caso contrario,
% elegir "no". Cualquier otra opción no imprime las graficas.
comparativa="si";


if apartado=="c" | apartado=="d"
    while t<0.002
    % Bucle correspondiente a la representacion del problema de Riemann
    % con solucion exacta.
        if apartado=="c"
            % Solucion del apartado c en la iteracion i-esima
            [U,P]=riemann(x,t,c0,rho,ur,ul,pr,pl,alf1R,alf2L);
        elseif apartado=="d"
            % Solucion del apartado d en la iteracion i-esima
            [U,P]=cauchy(x,t,c0,rho);
        end
        % Grafica correspondiente a las velocidades
        subplot(2,1,1);
        plot(x,U);
        %axis equal
        title("u(x,t) vs x")
        xlabel("x (m)")
        ylabel("u (m*s^{-1})")
        grid
        % Grafica correspondiente a las presiones
        subplot(2,1,2);
        plot(x,P);
        %axis([xmin xmax yminp ymaxp])
        title("p(x,t) vs x")
        xlabel("x (m)")
        ylabel(" p (N*m^{-2})")
        grid
        t=t+Dt;
        pause(Dt)
    end
elseif apartado=="er" | apartado=="ec"      % Ejecuta esta parte de codigo 
    u_i_1=zeros(size(x));                    % en el apartado e.
    p_i_1=zeros(size(x));
    if apartado=="er"
        % Apartado e aplicado al problema de riemann
        for i=1:N
            if x(i)<=0
                u_i_1(i)=ul;
                p_i_1(i)=pl;
            elseif x(i)>0
                u_i_1(i)=ur;
                p_i_1(i)=pr;
            end
        end
    elseif apartado=="ec"			    
        [u_i_1,p_i_1]=cauchy(x,0,c0,rho);		% Solucion de Cauchy 
    end
    u_i=u_i_1;
    p_i=p_i_1;
    while t<0.002           %El tiempo maximo elegido es suficiente para ver el comportamiento de las funciones
        for i=2:N-1   
            u_i_1(i)=u_i(i)*(1-Dt*c0/Dx)-(0.5*Dt/Dx)*((p_i(i+1)-p_i(i-1))/rho - c0*(u_i(i+1)+u_i(i-1)));
            p_i_1(i)=p_i(i)*(1-Dt*c0/Dx)-(0.5*Dt/Dx)*((u_i(i+1)-u_i(i-1))*rho*c0^2 -c0*(p_i(i+1)+p_i(i-1)));
        end
        u_i_1(1)=u_i_1(2);
        u_i_1(N)=u_i_1(N-1);
        p_i_1(1)=p_i_1(2);
        p_i_1(N)=p_i_1(N-1);
        if comparativa=="si"
            if apartado=="er"    
                [U,P]=riemann(x,t,c0,rho,ur,ul,pr,pl,alf1R,alf2L);
            elseif apartado=="ec"
                [U,P]=cauchy(x,t,c0,rho);
            end
        end
        subplot(2,1,1);
        if comparativa=="si"
            plot(x,u_i_1,x,U);
            legend("Aproximacion","Exacta")
        elseif comparativa=="no"
            plot(x,u_i_1);
        else 
            fprintf("No has seleccionado una opcion correcta")
            break
        end
        %axis equal
        title("u(x,t) vs x")
        xlabel("x (m)")
        ylabel("u (m*s^{-1})")
        grid
        % Grafica correspondiente a las presiones
        subplot(2,1,2);
        if comparativa=="si"
            plot(x,p_i_1,x,P);
            legend("Aproximacion","Exacta")
        elseif comparativa=="no"
            plot(x,p_i_1);
        else 
            fprintf("No has seleccionado una opcion correcta")
            break
        end
        %axis([xmin xmax yminp ymaxp])
        title("p(x,t) vs x")
        xlabel("x (m)")
        ylabel(" p (N*m^{-2})")
        grid
        t=t+Dt;
        pause(Dt)
        u_i=u_i_1;
        p_i=p_i_1;
    end
else 
    fprintf("No has seleccionado una opcion correcta")
end


% Funciones de calculo del problema
% Funcion del problema de Riemann

function [U,P]=riemann(x,t,c0,rho,ur,ul,pr,pl,alf1R,alf2L)
    %%% Funcion de velocidades del problema de Riemann para el sistema
    %%% hiperbolico de ondas.
    n=size(x,2);
    U=zeros(n,0);
    P=zeros(n,0);
    for i=1:n
        if x(i)+c0*t<=0
            U(i)=ul;
            P(i)=pl;
        elseif (x(i)+c0*t)*(x(i)-c0*t)<0
            U(i)=alf1R+alf2L;
            P(i)=(-alf1R+alf2L)*rho*c0;
        elseif x(i)-c0*t>0      
            U(i)=ur;
            P(i)=pr;
        end
    end 
end

% Funcion del problema de Cauchy que no es de Riemann para el apartado d.
function [U,P]=cauchy(x,t,c0,rho)
    % Pulsos correspondientes a los invariantes
    T=zeros(size(x));       % Array de tiempos necesario para las operaciones posteriores.
    T(:,:)=t
    v1_pulso=-exp(-(x+c0*T).^2);
    v2_pulso=exp(-(x-c0*T).^2);
    U=v1_pulso+v2_pulso;
    P=(-v1_pulso+v2_pulso)*rho*c0;
end

