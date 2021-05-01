% Media y covarianza de los brownianos

L=250;
t = linspace(0,5,L);

mu=1;
sig=0.6;
s=t(50);

M1=media_b_d(mu,t);
C1=cov_b_d(sig,t,s,L);

figure(1)
subplot(2,1,1)
plot(t,M1)
subplot(2,1,2)
plot(t,C1)

% Media y covarianza del browniano con deriva
function M=media_b_d(mu,t)
    M=mu*t;
end

function C=cov_b_d(sig,t,s,L)
    C=zeros(size(s,2));
    for i=1:L
        if t<s
            C(i)=sig^2*t;
        else
            C(i)=sig^2*s;
        end
    end
end

