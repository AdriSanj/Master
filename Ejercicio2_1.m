%%% Brownianos

L=250;
X = randn(1,L);
t = linspace(0,5,L);

mu=1;
sig=0.6;

Y1=brow_deriv(t,X,mu,sig,L);
Y2=brow_geo(t,X,mu,sig,L);
Y3=puen_brow(t,X,L);

figure(1)
plot(t,Y1)
title('Browniano con deriva')


figure(2)
plot(t,Y2)
title('Browniano geometrico')


figure(3)
plot(t,Y3)
title('Puente browniano')


% Browniano con deriva
function Y = brow_deriv(t,X,mu,sig,L)
    Y = zeros(1,L);
    for i = 1:L
        Y(i) = mu*t(i)+sig*X(i);
    end
end

% Browniano geometrico
function Y = brow_geo(t,X,mu,sig,L)
    Y = zeros(1,L);
    for i = 1:L
        Y(i) =exp(mu*t(i)+sig*X(i));
    end
end

% Puente browniano
function Y = puen_brow(t,X,L)
    Y = zeros(1,L);
    for i = 1:L
        Y(i) =X(i)-t(i)*X(1);
    end
end