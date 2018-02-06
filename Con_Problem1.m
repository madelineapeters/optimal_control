%This is the Script for the set up for Simple Conversion Rate Problem

global tau omega A K p B mu muz muG Pn
tau = 1; omega = 2; A = 370000; K = 370000*8500000/(370000-0.025*8500000); 
p = 0.000004; B = 16;
mu = 0.025; muz = 48; muG = 4; Pn = 43.859; 

%Hamiltonian
% H = G + LR*(A*(1 - (R/K)) - mu*R - p*R*M) + ...
% LM*(B*(1 - Ud)*p*Rd*Md*exp(-mu*tau) - muz*M - p*R*M + B*Pn*betapdf(t, 1,
% 1)*exp(-mu*t)) ...
% LG*(Ud*p*Rd*Md*exp(-mu*omega) - muG*G)
%The State, Adjoint, and Control functions needed for both processes
state1 = @(R,M) A*(1 - (R/K)) - mu*R - p*R*M;
state2 = @(t,R,Rd,M,Md,Ud) B*(1 - Ud)*p*Rd*Md*exp(-mu*tau) - muz*M - p*R*M + B*Pn*betapdf(t, 1, 1)*exp(-mu*t);
state3 = @(Rd,Md,G,Ud) Ud*p*Rd*Md*exp(-mu*omega) - muG*G;
adjoint1 = @(t,T,m1,m2,LR,LM,LMa,LGa,R,M,U) -LR*(-(A/K) - mu - p*M) - LM*(-p*M) ...
            - xalpha(t,T,tau)*LMa*B*(1 - U)*p*M*exp(-mu*tau) ....
            - xalpha(t,T,omega)*LGa*(U*p*M)*exp(-mu*omega);
adjoint2 = @(t,T,m1,m2,LR,LM,LMa,LGa,R,M,U) -LR*(-p*R) - LM*(- muz - p*R) ...
            - xalpha(t,T,tau)*LMa*B*(1 - U)*p*R*exp(-mu*tau) ...
            - xalpha(t,T,omega)*LGa*(U*p*R*exp(-mu*omega));
adjoint3 = @(LG) -1 - LG*(-muG);

%Function that choose how the new u compares to the old u
u_func = @(x,y)(x+y)/2;
