function [CD_L,cdf,Cd_Min,Cd_total,Total_Drag]=drag(tc,tc_lo,Tire_w,tire_D,y_w,...
    cord_W,high_wing,wing_span,Twist_angel,Alpha_induce_W,WSSL,density,...
    v_cruise,air_viscosity,i_w,Alpha_FOP,R2,RC_w,S_w,y_t,cord_T,R1,...
    high_tail,RC_Ht,RC_Vt,TC_Vt,S_Vt,length,lt,X_cg_dash,C_MGC,fuselage,...
    CL_wing_c,e_w,CL_tail,e_t,AR_w,AR_Ht,S_Ht,H_tail_span,alpha_twist,...
    TSSL,Alpha_induce_T,alpha_tail,Alpha_FOP_t)

%% friction drag

%(1)wing friction drag
% Mixed Boundary Layer 

[Cf_w,Swett_W] = Friction_Coff(y_w,wing_span,Twist_angel,Alpha_induce_W,WSSL,...
    density,cord_W,v_cruise,air_viscosity,i_w,Alpha_FOP,2*R2,high_wing,RC_w,1);
Cdf_w = (Swett_W/S_w)*Cf_w;

%(2) Horizontal Tail friction drag  
[Cf_T,Swett_T] = Friction_Coff(y_t,H_tail_span,alpha_twist,Alpha_induce_T,TSSL,...
    density,cord_T,v_cruise,air_viscosity,alpha_tail,Alpha_FOP_t,...
    R1*2,high_tail,RC_Ht,2);
Cdf_ht = (Swett_T/S_w)*Cf_T;

% (3)V tail friction drag
Re_vtr =  density*RC_Vt*v_cruise / air_viscosity;
Re_vtt =  density*TC_Vt*v_cruise / air_viscosity;
Cf_vtr = 1.328/sqrt(Re_vtr) ;
Cf_vtt = 1.328/sqrt(Re_vtt) ;
Cf_vt = (Cf_vtr+Cf_vtt)/2;
Svt_wet = 1.07*S_Vt*2;
Cdf_vt = Svt_wet/S_w * Cf_vt;
%% Interference factor 
Cdf_ht = Cdf_ht*1.05;
Cdf_vt = Cdf_vt *1.05;
%% Form Factor
FF_w = 1 + 2.7*(tc(1))+100*(tc(1))^4;                                          %suggested for the airfoil with t_c less than 21%
Cdf_w =  Cdf_w*FF_w;
FF_ht = 1 + 2.7*(tc(2))+ 100*(tc(2))^4;
Cdf_ht =Cdf_ht * FF_ht;
FF_vt = 1 + 2.7*(tc(3)) + 100*(tc(3))^4;
Cdf_vt = Cdf_vt*FF_vt;

%% fueslage Geometry
  [D_P,L_P,S_p,D_Cyl,L_Cyl,S_Cyl,D_CO,L_CO,S_CO,D_fL,D_fs,L_f,S_f]...
    = Basic_Geom(R2,R1,length,RC_w,X_cg_dash,C_MGC,lt);
                                           
A_WR_airfoil = 2*((tc_lo(1)+3)*(RC_w)^2*tc(1)/6) ;                          % wing root airfoil area*2 
A_HTR_airfoil = 2*((tc_lo(2)+3)*(RC_Ht)^2*tc(2)/6) ;                        % HT root airfoil area*2 
A_VTR_airfoil = ((tc_lo(3)+3)*(RC_Vt)^2*tc(3)/6) ;                          % VT root airfoil area 
if strcmp(fuselage,'jet') 
    % (1) jet transport fuselage(( parabloid nose + circular cylinderical middle portion + cone rea part))
    Sfues_wet = S_p + S_Cyl + S_CO-(A_WR_airfoil+A_HTR_airfoil+A_VTR_airfoil);
    S_pfus = 2/3*(D_P*L_P) +  D_Cyl*L_Cyl + 0.5*D_CO*L_CO;                  % fuselage planform area
    % (2) tadpole fuselage ( parabloid nose + circular cylinder + left frustum + right frustum)
elseif strcmp(fuselage,'tadpole')
    Sfues_wet = S_p + S_Cyl + sum(S_f)-(A_WR_airfoil+A_HTR_airfoil+A_VTR_airfoil);
    S_pfus = 2/3*(D_P*L_P) + D_Cyl*L_Cyl + sum((D_fL + D_fs).* L_f);        % fuselage planform area
end
%% fuelage friction drag
f = length/(2*R2);
Sf_max = (2*R1)^2 *pi/4;                                                    % maximum frontal area
FF_fues = 1+1.5/f^1.5+7/f^3;
Re_fues =  density*length*v_cruise / air_viscosity;
Cf_ffues = 0.455/(log10(Re_fues))^(2.58) ;
Cd_ffues_0e = FF_fues*Sfues_wet/S_w *Cf_ffues;                              % fuselage zero lift drag coff except base 
Cd_b_fuss = (0.029*(R1/R2)^3/(sqrt(Cd_ffues_0e*(S_w/Sf_max))))*(Sf_max/S_w); % fuselage zero lift drag coff base
Cd_ffues = Cd_ffues_0e + Cd_b_fuss;                                          % fuselage total zero-lift coff

%% drag coff due to lift
% (1) wing 
Cdi_w = CL_wing_c^2/(pi*e_w*AR_w);
% (2) tail
Cdi_t = (CL_tail^2/(pi*AR_Ht*e_t))*(S_Ht/S_w);
% (3) fuelsage 
Sb_fuss = (pi*(2*R1)^2)/4;                                                  % fuselage base area
C_dc = 1.2;                                                                 % experimental ss cross sectional drag of circular cylinder M<0.3
if   f >=2 && f <=4
    ef = 0.575;
elseif f >4 && f <=8
    ef = 0.625;
elseif f >8 && f <=12
    ef = 0.675;
elseif f >12 && f <=16
    ef = 0.725;
elseif f >16 && f <=20
    ef = 0.75;
elseif f >20 && f <=24
    ef = 0.775;
elseif f >24 && f <=28
    ef = 0.785;    
end
CD_fusL = 2*(-Alpha_FOP/57.3)^2 *(Sb_fuss/S_w) + ef*C_dc*...
    (abs(Alpha_FOP/57.3))^3 *(S_pfus/S_w);

%% Misc drag (tires and struts)
cd_tires = 3*Tire_w*tire_D/S_w*0.18;                                         % 3 because three tires in tricycle 
cd_strut = Tire_w*tire_D/S_w*1.2;
cdm = cd_tires + cd_strut;
%%
%CD_L = CD_fusL + Cdi_w + Cdi_t;
CD_L = Cdi_w + Cdi_t;
%cdf = Cdf_w + Cdf_ht + Cdf_vt+Cd_ffues;
cdf = Cdf_w + Cdf_ht + Cdf_vt;
%Cd_Min = (cdf+cdm)*1.25;
Cd_Min = (cdf)*1.25;
Cd_total = Cd_Min+CD_L;
Total_Drag = 0.5*density*(v_cruise^2)*S_w*Cd_total;
table(CD_L,cdf,Cd_Min,Cd_total,Total_Drag)
end