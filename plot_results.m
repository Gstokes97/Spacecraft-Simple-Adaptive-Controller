%Plot Results
%Graham Stokes &. 

%% q vs time plots
tt = t/60; 

figure(1)
hold on
subplot(3,1,1)
p1= plot(tt,rxyz.Data(:,1),'k',tt,rxyz.Data(:,2),'b',tt,rxyz.Data(:,3),'r')
xlabel('Time (Min)')
ylabel('MRP (SAC)')
legend('\sigma_1','\sigma_2','\sigma_3')
set(p1,'linewidth',0.35);
ylim([-0.5 1.2])
grid on
hold on

subplot(3,1,2)
p2=plot(tt,w_sc.Data(:,1),'k',tt,w_sc.Data(:,2),'b',tt,w_sc.Data(:,3),'r')
xlabel('Time (Min)')
ylabel('Angular Velocity (rad/s)')
legend('\omega_x','\omega_y','\omega_z')
ylim([-0.02 0.005])
set(p2,'linewidth',0.35);


grid on

subplot(3,1,3)
hold on
p3=plot(tt,tau.Data(:,1),'k',tt,tau.Data(:,2),'b',tt,tau.Data(:,3),'r')
xlabel('Time (Min)')
ylabel('Control Torque (Nm)')
legend('\tau_x','\tau_y','\tau_z')
ylim([-0.075 0.025])
set(p3,'linewidth',0.35);

grid on

