clear all
clc
global v0 CA0 k 


V_t0 = 0.5; % L
CA_t0 = 0.5; % mol/L
CB_t0 = 0.0; % mol/L
k = 8.4e-5; % 1/min
v0 = 0.002; % L/min
tF = 10*60; % min
CA0 = 0.05; %mol/L

tempo = [0 tF]; % min
CI = [V_t0 CA_t0 CB_t0];

[t,y] = ode23(@semibat, tempo, CI);

V  = y(:,1);
CA = y(:,2);
CB = y(:,3);

figure(1)
yyaxis left
plot(t,CA,'-bv','MarkerFaceColor','blue','MarkerSize',5,'LineWidth',2);
ylabel('C_A (mol/L)','FontSize',18)
hold on
yyaxis right
plot(t,CB,'--ro','MarkerFaceColor','red','MarkerSize',5,'LineWidth',2);
ylabel('C_B (mol/L)','FontSize',18)
xlabel('time (min)','FontSize',18)

hold off

figure(2)
plot(t,V,'--ro','MarkerFaceColor','red','MarkerSize',5,'LineWidth',2);
hold on
plot(t,V_t0+v0*t,'-b','MarkerFaceColor','blue','LineWidth',2);

ylabel('V (L)','FontSize',18)
xlabel('time (min)','FontSize',18)


function dydt = semibat(~,y)
    global v0 CA0 k 

    V  = y(1);
    CA = y(2);
    CB = y(3);
    
    neq = length(y);
    dydt = zeros(neq,1);
    
    dydt(1) = v0; % dV/dt
    dydt(2) = v0/V*CA0 - k*CA - CA/V*dydt(1); % dCA/dt
    dydt(3) = +2*k*CA - CB/V*dydt(1); % dCB/dt
    
end











