%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% heli.m
%% Matlab script to be run before Simulink files
%% ACS336 / ACS6336 / ACS6110
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;      % close all figures
clear all;      % clear workspace variables

%% Define Discrete Time MyDAQ Dynamics
T           = 0.015;            % Sample period (s)
ADC_Vres    = 20/((2^16)-1);    % ADC voltage resolution (V/bit)
Encoder_res = 2*pi/500;         % Encoder resolution (rad/wheel count)
DAC_Vres    = 20/((2^16)-1);    % DAC voltage resolution (V/bit)
DAC_lim_u   = 10;               % DAC upper saturation limit (V)
DAC_lim_l   = 0;                % DAC enforced lower saturation limit (V)
%% Define Continuous Time Helicopter Dynamics
g   = 9.81;     % Gravitational acceleration (ms^-2) 
% Rigid body parameters
% Masses and lengths
m1  = 0.0505;   % mass of fan assembly (kg)
m2  = 0.100;    % mass of counterweight (kg)
l1  = 0.110;    % distance from helicopter arm to elevation axis (m);
l2  = 0.070;    % distance from fan centres to pitch axis (m);
l3  = 0.108;    % distance from counterweight to elevation axis (m);
% Inertias
Je  = 2*m1*(l1^2)+m2*(l3^2);    % Inertia about elevation axis (kg*m^2);
Jt  = Je;                       % Travel axis inertia
Jp  = 2*m1*(l2^2);              % Pitch axis inertia
% Constraints
p_lim_u     = 80*pi/180;    % Upper pitch axis limit (rad)
p_lim_l     = -80*pi/180;   % Lower pitch axis limit (rad)
e_lim_u     = 50*pi/180;    % Upper elevation axis limit (rad)
e_lim_l     = -50*pi/180;   % Lower elevation axis limit (rad)

%% Ex 1: DETERMINE PITCH AXIS SPRING AND DAMPING COEFFICIENTS %%%%%%%%%%%%%
% Pitch axis spring and damping constants
k_s = 0.00038;           % Spring constant (kg*m^2*s^-2)
k_d = 0.00045;           % Viscous damping (kg*m^2*s^-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ex 2: DETERMINE POWER AMPLIFIER GAIN AND SATURATION LIMITS %%%%%%%%%%%%%
% Power amplifier
k_a         = 1.20;   % Power amplifier voltage gain
amp_sat_u   = 12;   % Power amplifier upper saturation limit (V)
amp_sat_l   = 0;   % Power amplifier lower saturation limit (V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Ex 3: CONSTRUCT FAN MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fan voltage - thrust steady state behaviour
V_ab   = 0:0.6:12;          % Fan voltage input (V)
Fss_ab = [0,0,0,0,0,0,0.0016,0.0029,0.0046,0.0064,0.0082,0.0107,0.0126,0.0155,0.0180,0.0206,0.0234,0.0261,0.0288,0.0315,0.0336];          % Steady-state fan thrust output (N)
% Fan voltage - thrust transient model.
FanDynData = load('FanDynData.txt');
TimeArray = FanDynData(1,:);
Excitation = FanDynData(2,:);
Excitation(Excitation==10) = Fss_ab(end);
Response = FanDynData(3,:);
% fan_dyn_data = iddata(Response',Excitation',T);
% fan_tf = tfest(fan_dyn_data,1);
opt = stepDataOptions('StepAmplitude',Fss_ab(end));
% step(fan_tf,TimeArray(1:sum(Excitation>0)),opt);
tau = 0.85;
fan_tf = tf(1,[tau 1]);
figure;
plot(TimeArray(1:sum(Excitation>0)),Response(Excitation>0),'Color',[0.3010 0.7450 0.9330])
hold on;
step(fan_tf,TimeArray(1:sum(Excitation>0)),opt,'b');
ylim([-0.005,0.055]);
%% Ex 4: DETERMINE EQUILIBRIUM CONTROL SIGNAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant control input required to achieve hover
U_e = 5.5;             % Simulated voltage output for hover
% U_e = 6.4;           % Voltage output from myDAQ

%% Ex 5: DEFINE LTI STATE-SPACE CONTROL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Approximate the fan voltage/thrust relationship by an affine           %
%  function of the form F_ab = alpha*V_ab+beta. Define alpha and beta.    %
p = polyfit(V_ab(7:end),Fss_ab(7:end),1);
alpha = p(1);
beta =  p(2);
figure;
plot(V_ab,Fss_ab,'kx');           % plot raw thrust data
grid on; hold on;
xlabel('Fan Voltage (V)');
ylabel('Output Thrust (N)');
plot(V_ab,alpha*V_ab+beta,'k-'); % plot linear approximation
%  State vector x:=[elev; pitch; trav; elev_dot; pitch_dot; trav_dot]     %
%  Note; these states model the dynamics of small perturbations around    %
%  the state of steady, level hover.                                      %
%  Define the control model given by x_dot = Ax + Bu, y = Cx + Du         %
A = [0 0                             0 1 0       0;
     0 0                             0 0 1       0;
     0 0                             0 0 0       1;
     0 0                             0 0 0       0;
     0 -k_s/Jp                       0 0 -k_d/Jp 0;
     0 -2*l1*(alpha*k_a*U_e+beta)/Jt 0 0 0       0];
B = [0               0;
     0               0;
     0               0;
     l1*alpha*k_a/Je l1*alpha*k_a/Je;
     l2*alpha*k_a/Jp -l2*alpha*k_a/Jp
     0               0];
C = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0];
D = [0 0;
     0 0;
     0 0];

%% Ex 6: Discrete Time Full state feedback control %%%%%%%%%%%%%%%%%%%%%%%%
% State feedback control design (integral control via state augmentation)
% x:=[elev; pitch; trav; elev_dot; pitch_dot; trav_dot; int_elev; int_trav]
% Define augmented system matrices
Cr      = [1 0 0 0 0 0
           0 0 1 0 0 0];    % Elevation and travel are controlled outputs
r       = 2;                                % number of reference inputs
n       = size(A,2);                        % number of states
q       = size(Cr,1);                       % number of controlled outputs
Dr      = zeros(q,2);
Aaug    = [A zeros(n,r); -Cr zeros(q,r)];
Baug    = [B; -Dr];
%Define LQR weighting matrices
% Qx = [2000 0 0 0 0 0 0 0;
%       0 20 0 0 0 0 0 0;
%       0 0 300 0 0 0 0 0;
%       0 0 0 200 0 0 0 0;
%       0 0 0 0 1.1 0 0 0;
%       0 0 0 0 0 100 0 0;
%       0 0 0 0 0 0 200 0;
%       0 0 0 0 0 0 0 0.3];    % State penalty
Qx = [1000 0 0 0 0 0 0 0;
      0 100 0 0 0 0 0 0;
      0 0 1 0 0 0 0 0;
      0 0 0 1000 0 0 0 0;
      0 0 0 0 10000 0 0 0;
      0 0 0 0 0 1000 0 0;
      0 0 0 0 0 0 5 0;
      0 0 0 0 0 0 0 10];
% Qx = [10000 0 0 0 0 0 0 0;
%       0 100 0 0 0 0 0 0;
%       0 0 1 0 0 0 0 0;
%       0 0 0 5000 0 0 0 0;
%       0 0 0 0 10000 0 0 0;
%       0 0 0 0 0 1000 0 0;
%       0 0 0 0 0 0 10 0;
%       0 0 0 0 0 0 0 10];
% Qx = [10000 0 0 0 0 0 0 0;
%       0 0.001 0 0 0 0 0 0;
%       0 0 1 0 0 0 0 0;
%       0 0 0 150000 0 0 0 0;
%       0 0 0 0 10000 0 0 0;
%       0 0 0 0 0 1000 0 0;
%       0 0 0 0 0 0 10000 0;
%       0 0 0 0 0 0 0 10];
Qx = [100 0 0 0 0 0 0 0;
      0 0.001 0 0 0 0 0 0;
      0 0 1 0 0 0 0 0;
      0 0 0 6000 0 0 0 0;
      0 0 0 0 10000 0 0 0;
      0 0 0 0 0 1000 0 0;
      0 0 0 0 0 0 100 0;
      0 0 0 0 0 0 0 200];
Qu = [10 0;
      0 10];    % Control penalty
% Discrete-Time LQR synthesis
Kdtaug  = lqrd(Aaug,Baug,Qx,Qu,T);      % DT state-feedback controller
Kdt     = Kdtaug(:,1:n); Kidt = -Kdtaug(:,n+1:end);
% Discrete-Time Kalman Filter Design
sysdt = c2d(ss(A,B,C,D),T,'zoh');     % Generate discrete-time system
Adt   = sysdt.a; Bdt = sysdt.b; Cdt = sysdt.c; Ddt = sysdt.d;
%  Kalman filter design; x_dot = A*x + B*u + G*w, y = C*x + D*u + H*w + v
Gdt     = 1e-1*eye(n);
Hdt     = zeros(size(C,1),size(Gdt,2)); % No process noise on measurements
Rw      = [1 0 0 0 0 0;
           0 1 0 0 0 0;
           0 0 1 0 0 0;
           0 0 0 1 0 0;
           0 0 0 0 1 0;
           0 0 0 0 0 1];   % Process noise covariance matrix
Rv      = 1e-05*eye(3);   % Measurement noise covariance matrix
sys4kf  = ss(Adt,[Bdt Gdt],Cdt,[Ddt Hdt],T);
[kdfilt, Ldt] = kalman(sys4kf,Rw,Rv);     % Kalman filter synthesis

%% Output Files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uncomment these lines when ready to implement your feedback controller  %
k_amp       = k_a;
Vfan        = V_ab;
Thrust      = Fss_ab;

KF_Initial  = [e_lim_l 0 0 0 0 0];
Kdi_Initial = [0 0];

csvwrite('aMatrix.txt',kdfilt.a)
csvwrite('bMatrix.txt',kdfilt.b)
csvwrite('cMatrix.txt',kdfilt.c(size(C,1)+1:end,:))
csvwrite('dMatrix.txt',kdfilt.d(size(C,1)+1:end,:))
csvwrite('crMatrix.txt',Cr)
csvwrite('KdMatrix.txt',Kdt)
csvwrite('KdiMatrix.txt',Kidt)
csvwrite('FanChar.txt',[Vfan,Thrust])
csvwrite('ModelParameters.txt',[T,U_e,l1,l2,l3,Jp,Jt,Je,m1,m2,g,...
    k_amp, amp_sat_u, amp_sat_l, DAC_lim_u, DAC_lim_l,...
    e_lim_l, e_lim_u, p_lim_l, p_lim_u])
csvwrite('KF_Initial.txt',KF_Initial)
csvwrite('Kdi_Initial.txt',Kdi_Initial)
