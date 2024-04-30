function [W_Total, W_Empty, Wing_Area, Power] = Weight_WingArea_Power(W_Payload, V_Stall_const, V_Max, V_Cruise, AR, rho)

%%%%%%%%%%%%%%%%%%% Initial Data %%%%%%%%%%%%%%%%%%%
W_Total = 30;        % Total Weight in kg (Initial Guess)
W_Batteries = 3*9.81;     % Newton (Two 22000 mAh capacity batteries)
W_B_to_W_T = W_Batteries / W_Total;
W_Empty = 20;        % kg (By estimation)
W_E_to_W_T = W_Empty / W_Total;

g = 9.81;             % gravity Acceleration
%W_Payload = 10*g;    % in Newten
A = 1;              % Const1 for Homebuilt UAV(Raymer's Table 3.1)
C = -0.1;           % Const2 for Homebuilt UAV(Raymer's Table 3.1)
%aa = -4.6*10^-5; 
%bb = 0.68;

%%%%%%%%%%%%%%%%%%% Weight Calculations %%%%%%%%%%%%%%%%%%%
syms x y
[W_Total,W_E_to_W_T] = vpasolve([x==(W_Payload)/(1-y),y==A*(x^C)],[x,y]);   %Sadraey's Eq.4.5 & Raymer's Table 3.1
%[W_Total,W_E_to_W_T] = solve([x==(W_Payload)/(1-y),y==aa*x+bb],[x,y]);   %Sadraey's Eq.4.5 & Eq.4.26
W_Total = round(double(W_Total),2);
W_E_to_W_T = double(W_E_to_W_T);
W_Empty = round(W_E_to_W_T * W_Total,2);

W_Empty = 18*g;
W_Total = W_Empty + 10*g;

fprintf("Total Weight = %f N = %f kg\n",W_Total,W_Total/g)
fprintf("Empty Weight = %f N = %f kg\n",W_Empty,W_Empty/g)



%%%%%%%%%%%%%%%%%%% Matching Plot Data %%%%%%%%%%%%%%%%%%%

% Stall Speed
rho_0 = 1.225;                             % Density at sea level in kg/m3
%V_Max = 18;                                % Max Velocity in m/s
V_Stall = linspace(0,V_Max);               % Stall Velocity in m/s
%V_Stall_const = 12.5;                        % Stall Velocity in m/s
Cl_max = 1.365;                             % Maximum Lift Coefficient

W_to_S = 0.5*rho_0*(V_Stall.^2)*Cl_max;     % Sadraey's Eq. 4.31
W_to_S_Stall = 0.5*rho_0*(V_Stall_const.^2)*Cl_max;


% Max Speed
e = 0.7;                 % Oswald span efficiency factor (0.7-0.85)
%AR = 5;                  % Aspect Ratio (4-7)
K = 1/(pi*e*AR);         % Induced Drag Factor
%rho = 1.1117;
sigma = rho/rho_0;       % Density at 1000 m altitude in kg/m3 %%%%%%%
eta_p = 0.7;             % Propeller Efficiency (0.7-0.85)
CD_0 = 0.035;            % Zero Lift Drag Coeff %%%%%%%


W_to_P_Vmax = (eta_p)./((0.5*rho_0*(V_Max^3)*CD_0*(1./W_to_S))+(((2*K)./(rho*sigma*V_Max))*(W_to_S)));    % Sadraey's Eq. 4.56

%CD_0 = (((2*P_Max*eta_p)/V_Max)-((4*K*(W^2))/(rho*sigma*(V_Max^2)*S)))/(rho_0*(V_Max^2)*S);  %%%%%%%


%Take-Off Run
V_TO = V_Stall*1.2;      % Take-Off Velocity (1.1-1.3)
%T_TO = 0.5*P_Max/V_TO;  % Estimated Engine Thrust for fixed-Pitch Propeller
S_TO = 80;               % Take-Off run requirement in meters %%%%%%
CL_Cruise = (W_Total)/(0.5*1.6*(V_Cruise^2)*rho);         % Cruise lift coefficient for subsonic aircraft 
D_CL_flapTO = 0;         %  take-off flap lift coefficient (0.3-0.8)
CL_TO = CL_Cruise+D_CL_flapTO;      % Take-Off Lift Coef
CD_0LG = 0.01;                      % Landing gear drag coefficient(0.006-0.012)
CD_0HLDTO = 0;                      % HLD (flap) drag coefficient at take-off (0.003-0.008)
CD_0TO =CD_0+CD_0LG+CD_0HLDTO;      % zero-lift drag coefficient at take-off
CD_TO = CD_0TO+K*(CL_TO^2);         % Drag coefficient at take-off
V_R = V_Stall*1.1;                  % Aircraft speed at rotation (1.1-1.2)
mu = 0.05;                          % Friction coefficient for runway surface (0.03-0.05)
CL_R = (2*W_to_S)./(rho*(V_R.^2));    % Lift coefficient at take-off rotation
CD_G = CD_TO-mu*CL_TO;

W_to_P_STO = ((1-exp(0.6.*rho.*g.*CD_G.*S_TO.*(1./W_to_S))).*(eta_p./V_TO))./(mu-(mu+(CD_G./CL_R)).*(exp(0.6.*rho.*g*CD_G.*S_TO.*(1./W_to_S))));     % Sadraey's Eq. 4.76


%Rate of Climb (ROC)
ROC = 1;                            % Rate of Climb at 1 m/s %%%%%%
eta_p_ROC = 0.55;                   % Propeller efficiency rated at (0.5-0.6)
S_wet_to_S_ref = (1.4*4)/1.4;         % Estimated Area wetted (exposed to air) from SolidWorks
K_LD = 9;                           % Constant for fixed gear prop aircrafts
Lift_to_Drag_Max = K_LD*sqrt(AR/S_wet_to_S_ref);      % Raymer's Eq. 3.12, sadraey's Typical (6-14) %%%%%

W_to_P_ROC = 1./((ROC./eta_p_ROC)+sqrt((2./(rho*sqrt((3*CD_0)./K)))*(W_to_S)*(1.155./(Lift_to_Drag_Max*eta_p_ROC))));     % Sadraey's Eq. 4.89


%Ceiling (Service Ceiling of 1000 m = 3281 ft)

W_to_P_Ceiling = sigma./((ROC./eta_p_ROC)+sqrt((2./(rho*sqrt((3*CD_0)./K)))*(W_to_S)*(1.155./(Lift_to_Drag_Max*eta_p_ROC))));     % Sadraey's Eq. 4.100


%Plot
plot(W_to_S,W_to_P_Vmax)
hold on
plot(W_to_S,W_to_P_STO)
plot(W_to_S,W_to_P_ROC)
plot(W_to_S,W_to_P_Ceiling)
xline(W_to_S_Stall,'k--');
legend('W/P_V_m_a_x','W/P_S_T_O','W/P_R_O_C','W/P_C_e_i_l_i_n_g','V_S_t_a_l_l Requirement')
ylim([0 0.8])



%%%%%%%%%%%%%%%%%%% Wing Area & Power %%%%%%%%%%%%%%%%%%%%%

for i = 1:length(W_to_S)
    if (W_to_S(i) >= W_to_S_Stall)
        n2=i;
        break;
    end
end
WingAreaLimit = 1.6;                              % Area Maximum limit
for i = 1:length(W_to_S)
    if (W_to_S(i) >= W_Total/WingAreaLimit)
        n1=i;
        break;
    end
end
l = zeros(4,1);
ll = zeros(n2,1);
for i = n1:n2
    l(1) = (eta_p)./((0.5*rho_0*(V_Max^3)*CD_0*(1./W_to_S(i)))+(((2*K)./(rho*sigma*V_Max))*(W_to_S(i))));
    l(2) = ((1-exp(0.6.*rho.*g.*CD_G.*S_TO.*(1./W_to_S(i)))).*(eta_p./V_TO(i)))./(mu-(mu+(CD_G./CL_R(i))).*(exp(0.6.*rho.*g*CD_G.*S_TO.*(1./W_to_S(i)))));
    l(3) = 1./((ROC./eta_p_ROC)+sqrt((2./(rho*sqrt((3*CD_0)./K)))*(W_to_S(i))*(1.155./(Lift_to_Drag_Max*eta_p_ROC))));
    l(4) = sigma./((ROC./eta_p_ROC)+sqrt((2./(rho*sqrt((3*CD_0)./K)))*(W_to_S(i))*(1.155./(Lift_to_Drag_Max*eta_p_ROC))));
    ll(i) = min(l);
end
[W_to_P_Final,index] = max(ll);
W_to_S_Final = W_to_S(index);

Wing_Area = W_Total/W_to_S_Final;
Power = W_Total/W_to_P_Final;

fprintf("Wing Area = %f m^2\n",Wing_Area)
fprintf("Power = %f watt\n",Power)

end
