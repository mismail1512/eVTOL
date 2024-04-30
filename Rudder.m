
% Rudder(40*0.514,1.5,28.5, 365,...
%     60,63,5.5,132*0.514,4.5,50,8,6.25,1.225,0.062 ,0,0.75,1.35,0,0,0.97)

function [Sr,Br,Cr]=Rudder(V_w,d_mostAft_mosFor,d_ACv_CgAf, Sw,...
    Bw,Lf,Df,Vstall,Clav,Sv,Bv,Cv,ruo,Vv,sidewash_gradient)
fprintf("***********Rudder **************\n");
V_w=V_w*2
%nv= 0.9;
%sidewash_gradient=0;
% Kf1=0.75 ;
% Kf2=1.35;
V_app=1.1*Vstall;
fprintf("The aircraft approach speed =%8.3f\n", V_app);
V_tot=sqrt(V_app^2+V_w^2);
fprintf("The aircraft total speed =%8.3f\n", V_tot);
Sr_Sv=0.25;
Sr=Sr_Sv*Sv;
fprintf("Rudder area =%8.3f m \n ", Sr);
Br_Bv=1;
Br=Br_Bv*Bv;
fprintf("Rudder span =%8.3f m \n ", Br);
CR_CV=0.3;
Cr=CR_CV*Cv;
fprintf("Rudder Chord =%8.3f m\n", Cr);
Ss=1.02*(Lf*Df+Sv);
%fprintf("The overall projected side area of the aircraft  =%8.3f m2\n", Ss)

Xca=((Lf*Df*0.5*Lf+Sv*((Lf-Cv)+0.5*Cv)))/(Lf*Df+Sv);
%fprintf("The overall center of the projected side area of the aircraft from the fuselage nose =%8.3f\n", Xca)
XAcVt_CgFor=d_ACv_CgAf-d_mostAft_mosFor;
Xcg=Lf-XAcVt_CgFor-0.75*Cv; % the distance between the center of gravity of the aircraft and the fuselage nose
dc=Xca-Xcg; %the distance between the center of the projected side area of the aircraft and the aircraft cg 
if (dc>0)
%fprintf("The center of the projected side area of the aircraft is  =%8.3f m behind the aircraft cg.\n", dc)
else
fprintf("The center of the projected side area of the aircraft is  =%8.3f m forward the aircraft cg.\n", dc)
end

Kf1=0.75;   
Kf2=1.35;
Cn0 =0;
Cy0 =0 ;
nv =0.97;
CDy =0.5;

Fw=0.5*ruo*Ss*CDy*V_w^2;  %CDy side drag coefficient--> 0.5–0.8.
%fprintf("the aircraft side force produced by the cross-wind =%8.3f N\n", Fw);
B=atan(V_w/V_app)
fprintf("the aircraft sideslip angle  =%8.3f Deg \n", B*57.3);
%Calculate the aircraft sideslip derivatives Cnβand Cyβ.


%CLα_v denotes the vertical tail lift curve slope
%dσ/dβ is the vertical tail sidewash gradient,
%ηv is the dynamic pressure ratio at the vertical tail
%The parameter Kf1 represents the contribution of the fuselage to the
%aircraft Cnβ  and depends strongly on the shape of the fuselage and its projected side area.
%Kf1 --> 0.65–0.85
%Kf2 --> 1.3–1.4
% The major contributor to the static directional stability derivative (Cnβ) --> A higher value for Cnβ implies a more directionally statically stable aircraft.
CnB=Kf1*Clav*(1-sidewash_gradient)*nv*XAcVt_CgFor*Sv/(Bw*Sw);
CyB=-Kf2*Clav*(1-sidewash_gradient)*nv*Sv/(Sw);

% Calculate the aircraft control derivatives CyδR and CnδR :
%the rudder angle of attack effectiveness (τ r) is extracted from Figure 12.12 for a control surface chord/lifting surface chord of 0.25.
tr=-6.624*(CR_CV^4)+12.07*(CR_CV^3)-8.292*(CR_CV^2)+3.295*(CR_CV)+0.004942;

Cydr=Clav*nv*tr*(Br*Sv)/(Sw*Bv);
%Vv the vertical tail volume coefficient
Cndr=-1*Clav*nv*Vv*tr*(Br)/(Bv);
syms def sidewash_angle
eqn1 = 0.5*ruo*(V_tot^2)*Sw*Bw*(Cn0+CnB*(B-sidewash_angle)+Cndr*def)+Fw*dc*cos(sidewash_angle)==0;
eqn2 = 0.5*ruo*(V_w^2)*Ss*CDy==0.5*ruo*(V_tot^2)*Sw*(Cy0+CyB*(B-sidewash_angle)+Cydr*def);
sol = solve([eqn1, eqn2], [def, sidewash_angle]);
deff = sol.def
sid = sol.sidewash_angle
defl=deff*180/pi;
SidewashAngle=sid*180/pi;
fprintf("The required rudder deflection (%8.3f deg) \n", defl);
fprintf("The sidewash angle (%8.3f deg) \n", SidewashAngle);

if (abs(defl)>30)
   fprintf("the rudder deflection is more than 30 deg, it is suggested to increase the rudder-tovertical tail chord ratio "); 
end
if (tr>1)
   fprintf(" If the angle of attack effectiveness of the rudder (τ R) is greater than 1 "); 
   fprintf(" there is no rudder which can satisfy the most critical directional control/trim requirement with the current vertical tail/aircraft cg combination.");
end


end

