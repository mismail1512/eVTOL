clear 
clc
close all

%%%%%%%%%%%%%%%%%Inputs%%%%%%%%%%%%%%%%%%%%

g = 9.81;                            % m/s2
W_Payload = 10*g;                    % Newton
V_Stall_const = 17;                  % m/s
V_Max = 30;                          % m/s
V_Cruise = 25;                       % m/s
AR = 6; 
lambda = 0.8;
[~, ~, ~, rho] = atmosisa(1000);                %atmospheric properties at cruise altitude 1000 m  at kg/m3
air_viscosity = 1.758*10^-5;
Wing_airfoil = '64(2)-415';
global R1
global R2
R1 = 0.04;                           % m
R2 = 0.21;                           % m
I_xx = 5;
Iyy = 4;
Phai_des = 30;                       % deg  Roll angle desired
V_w = 5;                             % m/s  Wind Speed
d_mostAft_mosFor = 0.01;             % m  cg tolerance
Lf = 2;                              % m  Fuselage Length
sidewash_gradient = 0;               %#######
Xmg = 0;
Xcg = 0.1;
Cmacwf = -0.08;
Zmg = 0;
[Zd,Zcg,ZT] = deal(0.3);
AoAws = 16;
AoAhse0 = 14;


%Initial Sizing

[W_Total, W_Empty, Wing_Area, Power] = Weight_WingArea_Power(W_Payload, V_Stall_const, V_Max, V_Cruise, AR, rho);


%Wing Design

[wing_span,RC_w,TC_w,Twist_angel,i_w,CL_wing_c,alpha_cruise,cl_alpha_3d_wf,MAC,Alpha_induce_W,WSSL,C_MGC,e_w]= ...
    Wing_Design(W_Total,W_Payload, Wing_Area, V_Cruise, AR, lambda);


%Tail Design

[l_opt, H_tail_span,  RC_Ht, TC_Ht,  V_tail_span,  RC_Vt, TC_Vt, Cl_h_req, CL_tail, alpha_tail,...
    S_Vt, S_Ht, cl_alpha_3d_Ht,C_MGC_Vt,V_VT,C_MGC_Ht,i_T,X_cg_F_dash,e_t,Alpha_induce_T,TSSL] = ...
    Taildesign(Wing_Area,AR,rho,W_Total, W_Empty, V_Cruise , lambda , alpha_cruise);


% Landing Gear 

[Gear_hight, z_nose_gear , z_main_gear] = designLandingGear( l_opt );


%Control Surface Design

% [t2,b_A,C_A, Aa,yi,y0] = Aileronn(Wing_Area,lambda,S_Ht,S_Vt,V_Stall_const,cl_alpha_3d_wf,...
%     I_xx,wing_span,RC_w,Phai_des,MAC);

[Sr,Br,Cr]=Rudder(V_w,d_mostAft_mosFor,l_opt, Wing_Area,...
    wing_span,Lf,0.4*R2,V_Stall_const,cl_alpha_3d_Ht,S_Vt,V_tail_span,C_MGC_Vt,rho,V_VT,sidewash_gradient);

% [elevator_span,elevator_chord,elevator_area]= Elevator(C_MGC_Ht,H_tail_span,S_Ht,rho,V_Stall_const,...
%     MAC,Xmg,Xcg,W_Total/g,Power/(V_Stall_const*1.2),(0.05*MAC)+Xcg,Iyy,...
%     l_opt-Xcg,cl_alpha_3d_Ht,Wing_Area,AR,V_Cruise,Cmacwf,Zd,Zmg,Zcg,i_w,i_T,cl_alpha_3d_wf,...
%     Power/(V_Max),ZT,AoAws,AoAhse0,V_Max,lambda);
% 
% 
% %Drag Estimation
% 
% tc = [0.15 0.12 0.12];                                                    % max thickness/cord ratio [wing HT VT]
% tc_lo = [0.349 0.3 0.3];                                                    % max thickness/cord ratio location [wing HT VT]
% %tires dimensions
% Tire_w = 0.06;                                                              % Tire width
% tire_D = 0.12;                                                               % Tire diameter
% y_w = [0  wing_span/2];                                                     % segments wing stattions locations
% cord_W = [RC_w  TC_w];                                                      % segments wing stations cord length
% high_wing =1;                                                               % 1 to use high wing and 0 for mid wing
% y_t = [0  H_tail_span/2];                                                   % segments wing stattions locations
% cord_T = [RC_Ht TC_Ht];                                                     % segments wing stations cord length
% high_tail =0;                                                    
% fuselage = 'tadpole';                                                       % choose tadpole or jet 
% Alpha_FOP = 0;
% length = 2.8;                                                               % Fuselage Length
% alpha_twist = 0.000001;
% Alpha_FOP_t = 0;
% 
% 
% [CD_L,cdf,Cd_Min,Cd_total,Total_Drag]=drag(tc,tc_lo,Tire_w,tire_D,y_w,...
%     cord_W,high_wing,wing_span,Twist_angel,Alpha_induce_W,WSSL,rho,...
%     V_Cruise,air_viscosity,i_w,Alpha_FOP,R2,RC_w,Wing_Area,y_t,cord_T,R1,...
%     high_tail,RC_Ht,RC_Vt,TC_Vt,S_Vt,length,l_opt,X_cg_F_dash,C_MGC,fuselage,...
%     CL_wing_c,e_w,CL_tail,e_t,AR,AR,S_Ht,H_tail_span,alpha_twist,...
%     TSSL,Alpha_induce_T,alpha_tail,Alpha_FOP_t);


%Propulsion


%Trim Condition

%[alpha_trim_deg, delta_e_deg] = trim(rho, V_Cruise, Wing_Area, W_Total,cl_alpha_3d_wf);
