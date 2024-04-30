function [wing_span,Root_Chord_w,Tip_Chord_w,Twist_angel,i_w,CL_wing_c,alpha_wing_root,Cl_alpha_3d_wf,MAC,Alpha_induce_W,WSSL,C_MGC,e_w] = Wing_Design(W_Total,W_Payload, Wing_Area, V_Cruise, AR, lambda)

% Inputs
S_w = Wing_Area;                                                            % Wing area
AR_w = AR;                                                                  % wing aspect ratio
Takeoff_Weight = W_Total;                                                   % Takeoff Gross weight
Landing_Weight = W_Total-W_Payload;                                         % Landing Gross weight
avergae_weight = (Takeoff_Weight+Landing_Weight)/2;                         % cruise average weight
[~, a_cruise, ~, density_cruise] = atmosisa(1000);                          % atmospheric properties at cruise altitude 1000 m
q = 0.5*density_cruise*V_Cruise^2*S_w;                                      % dynamic pressure @ cruise
% air_viscosity = 1.789*10^-5;                                              % air viscosity @ cruising altitiude
kinematic_viscosity_cruise = 1.581e-5;                                      % air
Twist_angel = -0.1;                                                         % wing geometric twist angle
Cl_c = (2*avergae_weight)/(density_cruise*V_Cruise^2 *S_w);               % Required wing Cruising Lift
global R1
global R2
% R1 = 0.04;                                                                % fuselage radius @ HT MGC
% R2 = 0.15;                                                                % fuselage radius @ cg
a = a_cruise;                                                               % sound speed at cruise
%% wing geometry
% equations are mainly from general aviation book

% Taper Wing
[wing_span,MAC,Root_Chord_w,Tip_Chord_w,C_MGC] = Geometry(AR_w,S_w,lambda);
%Re = (V_Cruise*MAC)/kinematic_viscosity_cruise;                             % Mean aerodynamic chord Reynold number

% 3D Lift Curve Slpoe
Cl_alpha_2d = 1.8*pi*(1+0.8*(15/100)) ;                                     % Wing airfoil Lift_curve slope(1/Degree) eq. 5.7 in Sadraey

% from 2d to 3d (sec. 9.5.4)
M = V_Cruise/a;
beta = sqrt(1-M^2);                                                         % Mach number parameter (Prandtl-Glauert)
kappa = Cl_alpha_2d/(2*pi);                                                 % two-dimensional lift curve slope (1/rad)
% the following equation can be got from Dr. Haitham first lecture notes, I guess
% SW_c2 = atan((Root_Chord_w/4-Tip_Chord_w/4)/(wing_span/2));               % Mostafa Kareem's derivation to mid cord sweep angle (it's valid only for zero quarter cord angle)
SW_c2 = 0;                                                                  % we assume no sweep angle Λ
Cl_alpha_3d_w = (AR_w*2*pi) / ( 2 + ...
    sqrt( (AR_w*beta/kappa)^2 * (1+(tan(SW_c2)^2)/beta^2) + 4) );           % Polhamus equation 9.72 (in 1/rad)
% account for fuselage effect
K_wf = 1 + 0.025*(2*R2/wing_span) - 0.25*(2*R2/wing_span)^2 ;
Cl_alpha_3d_wf = K_wf*Cl_alpha_3d_w;

% 3D Lift coff. @ alpha = 0 CL_0
Cl_0_2d = 0.3412;                                                           % 2D zero_lift coff
alpha_0_2d = -2.9643;                                                          % Airfoil Zero_lift angle (Cl = 0)
Cl_0_3d = alpha_0_2d * Cl_alpha_3d_wf;                                      % 3D zero_lift coff (eq. 9.67)

% don't know these equations source yet
alpha_3d_ZL_ET = -0.42;                                                     % change of wing angle of attacl relative to twist
alpha_0_3d = alpha_0_2d + alpha_3d_ZL_ET*Twist_angel;                       % Wing zero lift angle

%%  steps to determine a suitable wing AOI sec 9.3.6
% step 1: fuselage optimum AOA
alpha_F_opt = 0;                                                         % Fuselage Incidence angle. (-) is nose up
% step 2: cruise AOA
alpha_cruise = (1/Cl_alpha_3d_wf)*...                                       % cruisng angle of attack
    (avergae_weight/q) + alpha_0_3d;
% step 3: wing AOI
washout_correction_factor = ((1+2*lambda)/(3+3*lambda))*Twist_angel;        % twist correction factor (eq. 9.56)
% i_w = alpha_cruise + alpha_F_opt - washout_correction_factor;             % wing setting/incidence angle (deg)
i_w = 5;
alpha_wing_root = alpha_cruise-washout_correction_factor;                   % wing angle of attack @root
%% linear lifting line theory
Cl_alpha_3d = Cl_alpha_3d_wf;                                               % Wing airfoil lift curve slope (1/deg)
[CL_wing_c,e_w,Alpha_induce_W,WSSL] = Linear_Lifting_Line_Method(AR_w, Root_Chord_w,...
    wing_span, lambda, i_w, ...                                             % Linear Lifting Line Method
    Twist_angel,Cl_alpha_3d, alpha_0_2d,1, alpha_F_opt, density_cruise, V_Cruise);                                                         % Owstavald eff factor
%% results
table(wing_span, Root_Chord_w,Tip_Chord_w, MAC,Twist_angel,i_w, CL_wing_c, alpha_wing_root,Cl_c)
end