function x = del_runge_kutta_4(t,R,M,G,u,d1,d2,N,ode1,ode2,ode3)
% 4th order classic Runge-Kutta Method for solving

% state1 = @(R,M);
% state2 = @(t,R,Rd,M,Md,Ud);
% state3 = @(Rd,Md,G,Ud);

% constant step
h=1/N;
h2=h/2;

for i=1+d2:N+d2
    RK1 = ode1(R(:,i),M(:,i)); %c
    MK1 = ode2(t(:,i),R(:,i),R(:,i-d1),M(:,i),M(:,i-d1),u(:,i-d1)); %c
    GK1 = ode3(R(:,i-d2),M(:,i-d2),G(:,i),u(:,i-d2)); %c
    
    RK2 = ode1(R(:,i)+h2*RK1,M(:,i)+h2*MK1); %c
    MK2 = ode2(t(:,i)+h2,R(:,i)+h2*RK1,R(:,i-d1)+h2*RK1,M(:,i)+h2*MK1,M(:,i-d1)+h2*MK1,(u(:,i-d1)+u(:,i+1-d1))/2); %c
    GK2 = ode3(R(:,i-d2)+h2*RK1,M(:,i-d2)+h2*MK1,G(:,i)+h2*GK1,(u(:,i-d2)+u(:,i+1-d2))/2); %c
    
    RK3 = ode1(R(:,i)+h2*RK2,M(:,i)+h2*MK2); %c
    MK3 = ode2(t(:,i)+h2,R(:,i)+h2*RK2,R(:,i-d1)+h2*RK2,M(:,i)+h2*MK2,M(:,i-d1)+h2*MK2,(u(:,i-d1)+u(:,i+1-d1))/2); %c
    GK3 = ode3(R(:,i-d2)+h2*RK2,M(:,i-d2)+h2*MK2,G(:,i)+h2*GK2,(u(:,i-d2)+u(:,i+1-d2))/2); %c
    
    RK4 = ode1(R(:,i)+h*RK3,M(:,i)+h*MK3); %c
    MK4 = ode2(t(:,i)+h,R(:,i)+h*RK3,R(:,i-d1)+h*RK3,M(:,i)+h*MK3,M(:,i-d1)+h*MK3,u(:,i+1-d1)); %c
    GK4 = ode3(R(:,i-d2)+h*RK3,M(:,i-d2)+h*MK3,G(:,i)+h*GK3,u(:,i+1-d2)); %c
    
    R(i+1)=R(i)+(h/6)*(RK1+2*RK2+2*RK3+RK4); %c
    M(i+1)=M(i)+(h/6)*(MK1+2*MK2+2*MK3+MK4); %c
    G(i+1)=G(i)+(h/6)*(GK1+2*GK2+2*GK3+GK4); %c
end
x(1,:) = R;
x(2,:) = M;
x(3,:) = G;