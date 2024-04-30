function [l_opt, H_tail_span,  RC_Ht, TC_Ht,  V_tail_span,  RC_Vt, TC_Vt, Cl_h_req, CL_tail, alpha_tail, S_Vt, S_Ht,cl_alpha_3d_Ht,C_MGC_Vt,V_VT,C_MGC_Ht,i_T,X_cg_F_dash,e_t,Alpha_induce_T,TSSL] = Taildesign(S_w,AR_w,density,Total_Weight,Empty_Weight, v_cruise , lambda , alpha_wing)

global R1
global R2

%% wing 
% wing Airfoil    ACA 64(2)-415
Cl_Alpha_wing = 0.63335;
Cm_0_wing = -0.078 ; 
alpha_fueselage = 0;           %             
AR_w = 6;                                                                   % wing aspect ratio
sweep_angle = 0;      % review                                              % quarter cord sweep angle 
alpha_twist_wing = 0;

%% Aircraft Design Parameters and Conditions

%D_f = R1 + R2;
D_f = 2 * R2;

%% Assumption  .. 
X_cg_F_dash = 0.68;            % youssef         #######                           % Normalized forwar Centre of gravity(Refrence:wing LE)
X_ac_dash = 0.25;              % youssef         #######

%% wing geometry
[wing_span,C_avg,RC_w,TC_w,C_MGC] = Geometry(AR_w,S_w,lambda);

%% tail arm & Centre of gravity 
V_HT = 0.8;                          %review                                % Volume ratio for horizontal tail 
l_opt = 1.2 * ((4*S_w * ( V_HT* C_MGC))/(pi * (D_f)))^0.5 ;                 % Tail arm 

%% Lift Coefficient   
average_weight = (Total_Weight+Empty_Weight)/ 2; 
cl = (2*average_weight)/(density*v_cruise^2 *S_w);

%% Horizontal Tail angle of attack
T_eff = 0.96 ;                           % review                           % Horizontal Tail efficiency 
cm_0_w = ((AR_w*cos(sweep_angle)^2)/(AR_w + 2*cos(sweep_angle)))*(Cm_0_wing) + 0.01 * alpha_twist_wing ;
Cl_h_req = (cm_0_w+cl*(X_cg_F_dash-X_ac_dash))/(V_HT * T_eff);              % Required tail lift 

%% Tail Selected Airfoil 
% Must be symmetric and thinner than wing airfoil... NACA 0012 is selected
AR_Ht= (2/3) * AR_w;                                                           % Horizontal tail Aspect ratio

%% geometry of horizontal tail
S_Ht=(V_HT*S_w*C_MGC)/l_opt;                                                % horizontal tail area 

tr_Ht = 0.8;            % as wing
cl_alpha_2d_Ht = 1.8 * pi * (1 + 0.8 * 0.12);
cl_alpha_3d_Ht = cl_alpha_2d_Ht / (1 + (cl_alpha_2d_Ht / (pi * AR_Ht) )) ;

[H_tail_span,C_avg_Ht,RC_Ht,TC_Ht,C_MGC_Ht] = Geometry(AR_Ht,S_Ht,tr_Ht);

%% tail angle of attack in cruise 
alpha_tail = (Cl_h_req/cl_alpha_3d_Ht)* 180 / pi ;                          % Tail Angle of attack
alpha_tail_0 = alpha_tail;

%% the tail created lift coff
alpha_twist = 0.00001;                                                      % Tail Twist angle (deg)
alpha_0_TA = 0;                                                             % Tail Airfoil zero-lift angle of attack (deg)
Alpha_FOP_t = 0;

[CL_tail,e_t,Alpha_induce_T,TSSL] = Linear_Lifting_Line_Method(AR_Ht,RC_Ht,...
     H_tail_span,tr_Ht,alpha_tail,alpha_twist,cl_alpha_2d_Ht,alpha_0_TA,0,...
     Alpha_FOP_t,density,v_cruise);                                         % Linear Lifting Line Method

%% the Generated Lift Coefficient seems to be slightly less that the required " in magnitude ".Thus, adjusting the angle of attack is required
for i = 0 : 0.0001 : 2                                                     
    if (CL_tail - Cl_h_req)>0
        neg = -1;
    else
        neg = 1;
    end
[CL_tail,e_t,Alpha_induce_T,TSSL] = Linear_Lifting_Line_Method(AR_Ht,RC_Ht,...
     H_tail_span,tr_Ht,alpha_tail,alpha_twist,cl_alpha_2d_Ht,alpha_0_TA,0,...
     Alpha_FOP_t,density,v_cruise);                                         % Linear Lifting Line Method

alpha_tail = alpha_tail+i*neg;
  
 if (abs(CL_tail - Cl_h_req) <=0.001)
      break;
  end
end

%% downwash
epislon_0 = (2 * cl ) / (pi * AR_w );
epislon_alpha = (2 * Cl_Alpha_wing) / (pi * AR_w);
epislon = (epislon_0) + epislon_alpha*alpha_wing;

%% valculating Tail Incidence angle
i_T = alpha_tail - alpha_fueselage +epislon;                                    % horizontal tail incident angle

%% Displaying Data ..
% Horizontal Tail
disp( '______________________ Horizontal Tail ______________________  ' )
disp ('.');
disp ('   Tail arm    Span   Root_chord Tip_chord  Req_CL  Generat_CL  Alpha (deg)' )
disp ('===========================================================================')
disp ([l_opt H_tail_span  RC_Ht TC_Ht  Cl_h_req CL_tail alpha_tail ]) 

% checking for Cm_alpha 
Cm_alpha =   cl*(X_cg_F_dash-X_ac_dash) - cl_alpha_3d_Ht * T_eff * (S_Ht / S_w)*( (l_opt/C_avg_Ht) - X_cg_F_dash )*(1- epislon_alpha) ;
disp ('Checking the Contribution to the static longitudinal stability derivative (Cmα)');
fprintf ('(Must be negative to insure a stabilizing contribution)  ===>  Cm_alpha  =');
disp( Cm_alpha );


%% geometry of vertical tail
% conventional Configuration is selected
V_VT = 0.09;                                                                % volume ratio for vrtical tail 
AR_Vt=1.9;  
S_Vt=(V_VT*wing_span*S_w)/l_opt;                                            % vertical tail area
S_Vt = S_Vt/2;                                                              % as the vertic tail is only in one side 
tr_Vt = 0.8;                                                                % Taper Ratio of Vertical Tail. 
[V_tail_span,C_avg_Vt,RC_Vt,TC_Vt,C_MGC_Vt] = Geometry(AR_Vt,S_Vt,tr_Vt);


%% Vertical Tail Data
disp( '_______________ Vertical Tail _________________  ' )
disp ('.');
disp ('   Tail arm    Span   Root_chord Tip_chord  ' )
disp ('==============================================')
disp ([l_opt V_tail_span  RC_Vt TC_Vt ])


end   
