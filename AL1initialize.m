% QuadrotorL1initialize.m
%% define quadrotor constants
% mass, inertias & geometry (eqs 14 & 15)
m=24.55;     % mass (kg)
Ixx=1.224;  % I_x (kgm^s)
Iyy=6.688;  % I_y (kgm^s)
Izz=6.688; % I_z (kgm^s)
g=9.81;      % gravitational constant (ms^-2)
r=0.13;    % arm length (m)
d=0.57;    % arm length (m)
% motors & rotors (eqs 16 & 17)
Vmax=12.5;  % max voltage
kf=0.0031;    % thrust constant f = k_f*v^2 (eq 16)\approx g/90 

%% linearized open loop system (eqs 45 - 47)
% state variables = \delta  [U ,V ,W ,P ,Q ,R ,X ,Y ,Z ,\phi ,\theta , \psi]
A=  [0 0 0     0 0 0    0 0 0    0 -g 0;...
    0 0 0     0 0 0    0 0 0    g 0 0;...
    0 0 0     0 0 0    0 0 0    0 0 0;...
    0 0 0     0 0 0    0 0 0    0 0 0;...
    0 0 0     0 0 0    0 0 0    0 0 0;...
    0 0 0     0 0 0    0 0 0    0 0 0;...
    1 0 0     0 0 0    0 0 0    0 0 0;...
    0 1 0     0 0 0    0 0 0    0 0 0;...
    0 0 1     0 0 0    0 0 0    0 0 0;...
    0 0 0     1 0 0    0 0 0    0 0 0;...
    0 0 0     0 1 0    0 0 0    0 0 0;...
    0 0 0     0 0 1    0 0 0    0 0 0];
b=[  0     0     0     0;...
    0     0     0     0;...
    -1/m    0     0     0;...
    0    1/Ixx  0     0;...
    0     0   1/Iyy   0;...
    0     0     0 1/Izz;...
    0     0     0     0;...
    0     0     0     0;...
    0     0     0     0;...
    0     0     0     0;...
    0     0     0     0;...
    0     0     0     0];
C = [0 0 0     0 0 0    1 0 0    0 0 0;...
    0 0 0     0 0 0    0 1 0    0 0 0;...
    0 0 0     0 0 0    0 0 1    0 0 0;...
    0 0 0     0 0 0    0 0 0    0 0 1];
D=zeros(4,4);
C_full=[...
    0 0 0     0 0 0    1 0 0    0 0 0;...
    0 0 0     0 0 0    0 1 0    0 0 0;...
    0 0 0     0 0 0    0 0 1    0 0 0;...
    0 0 0     0 0 0    0 0 0    1 0 0;...
    0 0 0     0 0 0    0 0 0    0 1 0;...
    0 0 0     0 0 0    0 0 0    0 0 1];
x_initial=[0;0;0;0;0;0;0;0;0;0;0;0];

%% control allocation matrix (eq 19)
Lambda= [
    -kf         kf        -kf         kf;...
    -r*kf      -r*kf     -r*kf       -r*kf;...
    0         -d*kf   0          d*kf;...
    d*kf    0        -d*kf    0];
B=b*Lambda;
v_trim=m*g/4/kf;                          % trim voltage (eq 49)
V_trim=[v_trim; v_trim; v_trim; v_trim];  % trim voltage vector

%% LQR inner loop controller design - altitude and attitude
% stability and acontrollability check
sys_open=ss(A,B,C,D);
co=ctrb(sys_open);    % 6poles at origin
rank(co);             % rank =12, the system is controllable

%weighting matrices
Q=100 * eye(12) ;
R=eye(4);

K1 = lqr(A,B,Q,R);  %feed back gain (eq 50)

%% desired close loop system (Section 3.1)
Am = (A-B*K1);
Bm = B;
%% State predictor model (Section 3.2)
Bum=[...  % Bm' Bum = 0
    1     0     0     0     0     0     0     0;...
    0     1     0     0     0     0     0     0;...
    0     0     0     0     0     0     0     0;...
    0     0     0     0     0     0     0     0;...
    0     0     0     0     0     0     0     0;...
    0     0     0     0     0     0     0     0;...
    0     0     1     0     0     0     0     0;...
    0     0     0     1     0     0     0     0;...
    0     0     0     0     1     0     0     0;...
    0     0     0     0     0     1     0     0;...
    0     0     0     0     0     0     1     0;...
    0     0     0     0     0     0     0     1];
Cm = eye(12);
Dm=zeros(12,4);
% sys_cl=ss(Am,Bm,C,D) % closed inner loop system

%% Adaption laws (eqs 31-35 )
Gamma=1000;              % adaption gain
P=lyap(Am',eye(12));     % Lyapunov matrix
PBmG=P*Bm*Gamma;         % P B_m \Gamma
PBumG=P*Bum*Gamma;       % P B_um \Gamma

%% Controller laws (eqs 36-41)
Kg=-(C*Am^-1*Bm)^-1;     % tracking pre-compensator K_g
%HD=ss(eye(4)*tf([1],[1 0]))                % filter H_D(s) (tf form)
HD=ss(zeros(4,4),eye(4),eye(4),zeros(4,4)); % filter H_D(s) (ss form)
K2=eye(4)*10; % feedback gain matrix

Hxm=ss(Am,Bm,eye(12),zeros(12,4));    % eq 36
Hxum=ss(Am,Bum,eye(12),zeros(12,8));  % eq 36
Hm=ss(Am,Bm,C,D);                     % eq 37
Hum=ss(Am,Bum,C,zeros(4,8));          % eq 37
HmHum=(Hm^-1)*Hum;                    %  H_m^{-1} H_{um} (eq 40)
% get minimal and proper realization of H_m^{-1} H_{um}
% pole(HmHum), zero(HmHum) % poles and zeros
xx=zpk(HmHum);
for i=1:4
    for j=1:8
        if abs(xx(i,j).k)<1e-12,
            xx(i,j)=0;
        else
            xx(i,j)=minreal(xx(i,j));
        end
    end
end
% make proper
xxx2=zpk([],[-100, -200],1);
xxx3=zpk([],[-100, -200,-400],1);
xxx1=zpk([],-100,1);
HmHumprop=xx*append(xxx2,xxx2,xxx3,xxx3,xxx1,xxx1,xxx1,xxx1);
HmHum=ss(HmHumprop); %  H_m^{-1} H_{um}

%% Low pass filter design
HC=ss(-K2,eye(4),K2,zeros(4,4));

%% define failure level
vfail=(12)^(1/3.88) 

%% Projection bounds
theta1hatbound=40;
theta2hatbound=2;
sigama1hatbound=10;
sigama2hatbound=1;
Deltabound=1;
   
%% check L_1 bounds
if(0)
    L1HmHCKg=L1norm(Hm*HC*Kg)
    L1Fm=L1norm(Hxm*ss(-K2,eye(4),-K2,eye(4))) % || Hxm*(I-H_C) ||
    HxmHCinvHmC= Hxm*HC*(Hm^-1)*C;
    Fum =(eye(12)-HxmHCinvHmC)*Hxum;
    
end


