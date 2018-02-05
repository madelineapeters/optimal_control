function z = del_back_runge_kutta_4(t,T,LR,LM,LG,R,M,u,d1,d2,N,ode1,ode2,ode3)
% 4th order classic Runge-Kutta Method for solving
% like regular RK, but solved backward in time

%adjoint1 = @(t,T,m1,m2,LR,LM,LMa,LGa,R,M,U)
%adjoint2 = @(t,T,m1,m2,LR,LM,LMa,LGa,R,M,U)
%adjoint3 = @(LG)

% constant step
h=1/N;
h2=h/2;

for j=1:N
    i = N+d2+2-j;    
    
    LR1 = ode1(t(:,i),T,d1,d2,LR(:,i),LM(:,i),LM(:,i+d1),LG(:,i+d2),R(:,i),M(:,i),u(:,i));
    LM1 = ode2(t(:,i),T,d1,d2,LR(:,i),LM(:,i),LM(:,i+d1),LG(:,i+d2),R(:,i),M(:,i),u(:,i));
    LG1 = ode3(LG(:,i));
    
    LR2 = ode1(t(:,i),T,d1,d2,LR(:,i)-h2*LR1,LM(:,i)-h2*LM1,LM(:,i+d1)-h2*LM1,LG(:,i+d2)-h2*LG1,(R(:,i)+R(:,i-1))/2,(M(:,i)+M(:,i-1))/2,(u(:,i)+u(:,i-1))/2);
    LM2 = ode2(t(:,i),T,d1,d2,LR(:,i)-h2*LR1,LM(:,i)-h2*LM1,LM(:,i+d1)-h2*LM1,LG(:,i+d2)-h2*LG1,(R(:,i)+R(:,i-1))/2,(M(:,i)+M(:,i-1))/2,(u(:,i)+u(:,i-1))/2);
    LG2 = ode3(LG(:,i)-h2*LG1);
    
    LR3 = ode1(t(:,i),T,d1,d2,LR(:,i)-h2*LR2,LM(:,i)-h2*LM2,LM(:,i+d1)-h2*LM2,LG(:,i+d2)-h2*LG2,(R(:,i)+R(:,i-1))/2,(M(:,i)+M(:,i-1))/2,(u(:,i)+u(:,i-1))/2);
    LM3 = ode2(t(:,i),T,d1,d2,LR(:,i)-h2*LR2,LM(:,i)-h2*LM2,LM(:,i+d1)-h2*LM2,LG(:,i+d2)-h2*LG2,(R(:,i)+R(:,i-1))/2,(M(:,i)+M(:,i-1))/2,(u(:,i)+u(:,i-1))/2);
    LG3 = ode3(LG(:,i)-h2*LG2);
    
    LR4 = ode1(t(:,i),T,d1,d2,LR(:,i)-h*LR3,LM(:,i)-h*LM3,LR(:,i+d1)-h*LM3,LG(:,i+d2)-h*LG3,R(:,i-1),M(:,i-1),u(:,i-1));
    LM4 = ode2(t(:,i),T,d1,d2,LR(:,i)-h*LR3,LM(:,i)-h*LM3,LR(:,i+d1)-h*LM3,LG(:,i+d2)-h*LG3,R(:,i-1),M(:,i-1),u(:,i-1));
    LG4 = ode3(LG(:,i)-h*LG3);
    
    LR(i-1)=LR(i)-(h/6)*(LR1+2*LR2+2*LR3+LR4); 
    LM(i-1)=LM(i)-(h/6)*(LM1+2*LM2+2*LM3+LM4); 
    LG(i-1)=LG(i)-(h/6)*(LG1+2*LG2+2*LG3+LG4); 
end
z(1,:) = LR;
z(2,:) = LM;
z(3,:) = LG;