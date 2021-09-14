%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPACECRAFT ATTITUDE DYNAMICS AND CONTROL - INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steve Ulrich
% 10 October, 2020
% AERO 4540
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all

disp('Initalizing Parameters...')

s = 0.01; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RW time constant 
rw_t = 0.01; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reaction Wheel 1 about Spin Axis a1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tex_B = [0 0 0]';
%Ja1        = 4.10e-4; % wheel 1 moment of inertia about its spin axis, kg m^2
TCt = [0 0 0]';  

Ja1 = 0.0106;
Ja1_inv    = 1/Ja1; % inverse of wheel 1 inertia, kg^-1 m^-2
wa1_ini    = 0; % wheel 1 initial angular rate about its spin axis wrt ECI, rad/s
ha1_ini    = Ja1*wa1_ini; % wheel 1 initial angular momentum about its spin axis, Nms
a1_B       = [1,0,0]'; % components of wheel 1 spin axis in BOF
ha1_B_ini  = ha1_ini*a1_B; % components of initial angular momentum vector of wheel 1 in BOF, Nms
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reaction Wheel 2 about Spin Axis a2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ja2        = 4.11e-4; % wheel 1 moment of inertia about its spin axis, kg m^2
Ja2 = 0.0106;
Ja2_inv    = 1/Ja2; % inverse of wheel 2 inertia, kg^-1 m^-2
wa2_ini    = 0; % wheel 2 initial angular rate about its spin axis wrt ECI, rad/s
ha2_ini    = Ja2*wa2_ini; % wheel 2 initial angular momentum about its spin axis, Nms
a2_B       = [0,1,0]'; % components of wheel 2 spin axis in BOF
ha2_B_ini  = ha2_ini*a2_B; % components of initial angular momentum vector of wheel 2 in BOF, Nms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reaction Wheel 3 about Spin Axis a3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ja3        = 4.12e-4; % wheel 1 moment of inertia about its spin axis, kg m^2
Ja3 = 0.0106;
Ja3_inv    = 1/Ja3; % inverse of wheel 3 inertia, kg^-1 m^-2
wa3_ini    = 0; % wheel 3 initial angular rate about its spin axis wrt ECI, rad/s
ha3_ini    = Ja3*wa3_ini; % wheel 3 initial angular momentum about its spin axis, Nms
a3_B       = [0,0,1]'; % components of wheel 3 spin axis in BOF
ha3_B_ini  = ha3_ini*a3_B; % components of initial angular momentum vector of wheel 3 in BOF, Nms

%% Model Initialization + constant

wnn = 0.02; 
zetaa = 1;
alpha = 1;

A = [zeros(3,3), eye(3); (-1*(wnn)^2)*eye(3),-2*zetaa*wnn*eye(3)];
B = [zeros(3,3);(wnn^2)*eye(3)];
C = [eye(6)];
D = [zeros(6,3)];


%% SAC Paper Inertia %%
gIE = 1e3*(eye(3)); 
gPE = 1e3*(eye(3));
gIX = 1e3*(eye(6));
gPX = 1e3*(eye(6));
gIU = 1e3*(eye(3));
gPU = 1e3*(eye(3));
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spacecraft Platform (without the wheels) Dynamics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jsc = 1*[114 0 0;0 86 0;0 0 87];

JscC = Jsc/1.5;
Jsc_inv = Jsc^-1;

wsc_ini_B   = [0,0,0]'; % components of initial spacecraft angular velocity vector wrt ECI in BOF (rad/s)
hsc_ini_B   = Jsc*wsc_ini_B; % components of initial spacecraft angular momentum vector in BOF (Nms)

%%%%%%%%%%%%%%%%%%
%% ATTITUDE INI %%
%%%%%%%%%%%%%%%%%%

%%%%%%%%%% For Quaternion %%%%%%%%%
a1 = 1;
a2 = 5;
a3 = -2;
aNorm = sqrt(a1^2 + a2^2 + a3^2); 

a = [a1/aNorm a2/aNorm a3/aNorm]; 

phi = 20; 
d2r = pi()/180; 
phiR = phi*d2r; 

qsc_I_ini = [a(1)*sin(phiR/2) a(2)*sin(phiR/2) a(3)*sin(phiR/2) cos(phiR/2)];

%%%%%%%%%%%%%% For MRP %%%%%%%%%%%%%
r_I_ini_pos= (qsc_I_ini(1:3)/qsc_I_ini(4));
%r_I_ini = [-0.1 0.5 1.0]; %Ulrichs Initial MRP
r_I_ini_dot = [0 0 0];
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Total Angular Momentum (spacecraft platform + wheels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

htot_B_ini  = hsc_ini_B + ha1_B_ini+ ha2_B_ini + ha3_B_ini; % components of initial total angular momentum vector in BOF, Nms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Attitude Sensors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wrel_noise = 2*pi/60;           % tachometer noise in rad/s (1 RPM)
wrel_var = wrel_noise^2;        % tachometer variance

qsc_noise = 1/(60*60) * pi/180; % star tracker noise in rad (1 arcsec)
qsc_var = qsc_noise^2;          % star tracker variance

wsc_noise = pi/180;             % rate gyro noise in rad/s (1 deg/s)
wsc_var = wsc_noise^2;          % rate gyro variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Control Gains and Desired Quaternion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphaP = 1; 
alphaM = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Desired Quaternion %%

q_d = [0 0 0 1]'; % desired quaternion




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Simulation and Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Running the sim with different gains
disp('Opening Simulation')
run_SACRW
disp('Simulation Complete')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('plotting results')
plot_results
sumtorques 





