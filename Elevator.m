function [elevator_span,elevator_chord,elevator_area]= Elevator(Ch,Bh,Sh,ruo,V_stall,...
    mean_chord,Xmg,Xcg,mass,thrust,Xacwf,Iyy,...
    lm,Clah,Sw,AR,Vc,Cmacwf,Zd,Zmg,Zcg,AoA_W,ih,Clawf,Tmax,ZT,AoAws,AoAhse0,Vmax,lambda)
%CLαh is the tail lift curve slope 
fprintf("***********Elevator **************\n");

%CE_CH=0.37;
BE_BH=1;
%Cmacwf=-0.08;
elevator_span=BE_BH*Bh;


deltaCl_flap=0;

%%%%%%%%%%%%%%%%%%%%%%%%%   Take off   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
max_Def=20 ; %deg
ang_acc= 12 ; %deg*s2
Oswald = 1.78*(1 - 0.045*AR^0.68) - 0.64; % e is the Oswald span efficiency factor
k=1/(3.14*Oswald*AR);    %where K is the wing induced drag factor
Clcruise=2*mass*9.81/(ruo* (Vc^2)*Sw);
%deltaCl_flap=0.5;
ClTO=Clcruise+deltaCl_flap;
%CdoTO = CDo + CDoLG + CDoHLD_TO
%where CDo is the clean aircraft zero-lift drag coefficient (see Table 4.12)--> 0.02–0.03,
CDo=0.025;
%CDoLG is the landing gear drag coefficient --> 0.006 to 0.012
CDoLG=0.008;
%CDoHLD_TO is the high-lift device drag coefficient at take-off configuration-->0.003 to 0.008
CDoHLD_TO=0;
CdoTO = CDo + CDoLG + CDoHLD_TO;
Cdtot=CdoTO+k*ClTO^2; %where the aircraft zero-lift drag coefficient at take-off configuration (CDoTO)
VR= V_stall; % Rotataion Speed  
DTO=0.5*ruo* (VR^2)* Cdtot*Sw; %aerodynamic drag during take off
LwfTO=0.5*ruo* (VR^2)* ClTO*Sw; %wing/fuselage lift  during take off
Macwf=0.5*ruo* (VR^2)* Cmacwf*Sw*mean_chord; %wing/fuselage pitching moment Cmacwf=0.05 ########
W= mass*9.81;    %Weight
fprintf("Weight =%8.3f N \n", W);

Ff=0.04*(W-LwfTO);  % friction force
acc=(thrust-DTO-Ff)/mass;
fprintf("Aircraft linear acceleration at the time of take-off rotation =%8.3f m/s^2 \n", acc);

Mw=abs(W*(Xmg-Xcg))  %aircraft weight moment  
Md=abs(DTO*(Zd-Zmg))  %aircraft drag momen  
Mt=abs(thrust*(ZT-Zmg))  %engine thrust moment  
MLwf=abs(LwfTO*(Xmg-Xacwf)) %wing/fuselage lift moment (MLwf),  
Ma=abs(mass*acc*(Zcg-Zmg))  % linear acceleration moment (Ma) ????

Lh=(MLwf+Macwf+Ma-Mw+Md-Mt-Iyy*ang_acc/57.3)/(lm);
fprintf("the desired horizontal tail lift (Lh) during take-off rotation =%8.3f\n", Lh);
% 
% Mlh=-Lh*(lm);   %horizontal tail lift moment (MLh)  ???    
% fprintf("the desired horizontal tail lift moment =%8.3f\n", Mlh);

Clh=(2*Lh)/(ruo*Sh*VR^2);   %the tail lift coefficient
fprintf("the desired horizontal tail lift coefficient =%8.3f\n", Clh);


 e0=2*Clcruise/(3.14*AR) %e0 is downwash angle at zero angle of attack % rad
 de_da=2*Clawf/(3.14*AR)   %deg/deg
 e=e0+de_da*AoA_W/57.3      %#####################
 %AoA_h=AoA_W+ih-e*57.3   %deg
 AoA_h=-7
n1=AoA_W
n2=ih
fprintf("the horizontal tail angle of attack at the instance of take-off rotation =%8.3f deg \n", AoA_h);

te=(((AoA_h/57.3)+(Clh/Clah))/(-1*max_Def/57.3));
fprintf("The angle of attack effectiveness of the elevator =%8.3f\n", te);

syms CE_CH
eqn = te==-6.624*(CE_CH^4)+12.07*(CE_CH^3)-8.292*(CE_CH^2)+3.295*(CE_CH)+0.004942;
sol = vpasolve([eqn],[CE_CH]);
if (isreal(sol(1)));
    A=sol(1);
else A=10  ; 
end
if (isreal(sol(2)));
    B=sol(2);
else B=10  ; 
end
if (isreal(sol(3)))
    C=sol(3);
else C=10   ;
end
if (isreal(sol(4)));
    D=sol(4);
else D=10   ;
end
r=[A B C D]

CE_CH = min(r)
elevator_chord=CE_CH*Ch;
elevator_area=elevator_chord*elevator_span;
fprintf("elevator area =%8.3f m2\n", elevator_area);
fprintf("elevator chord =%8.3f m\n", elevator_chord);
fprintf("elevator span =%8.3f m\n", elevator_span);
SE_SH=elevator_area/Sh;

% Te=-6.624*(CE_CH^4)+12.07*(CE_CH^3)-8.292*(CE_CH^2)+3.295*(CE_CH)+0.004942;
% fprintf("The effective angle of attack from the mathematical model  =%8.3f\n", Te);

% the elevator deflection (δE) required to maintain longitudinal trim at various flight conditions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eh=0.95;
lh=lm+Xcg;
VH=lh*Sh/(Sw*mean_chord);
CmdE=-Clah*eh*VH*BE_BH*te;
CldE=Clah*eh*SE_SH*BE_BH*te;

ClhdE=Clah*te;
Q=0.5*ruo*Vc^2;  %pa
xMgAc=0.4;
xMgCgAf= 0.2;
xActCgAf=1;
Cma=Clawf*(xMgAc-xMgCgAf)-Clah*eh*SE_SH*(xActCgAf/mean_chord)*(1- de_da);
Cl1=2*W/(ruo*Sw*Vc^2);
Cm0=  0.05;%######################
Cl0=  0.3412;%###################### 

dE=-(((Tmax*0/(Q*Sw*mean_chord))+Cm0)*Clawf+(Cl1-Cl0)*Cma)/(Clawf*CmdE-Cma*CldE);
fprintf(" the elevator deflection (δE) required to maintain longitudinal trim  =%8.3f  rad =%8.3f deg \n", dE,dE*57.3);

%the wing stall angle:
AoAwsTO=AoAws-AoA_W;
fprintf(" the wing stall angle: =%8.3f deg \n", AoAwsTO);

    
%the horizontal tail take-off angle
AoAhTO=AoAwsTO*(1-de_da)+ih-e0*57.3;
fprintf(" the horizontal tail take-off angle =%8.3f  deg \n", AoAhTO);

%The tail stall angle of attack during take-off rotation (αhs) 
%dAOAhE the magnitude of reduction in tail stall angle of attack due to elevator deflection and is determined using Table 12.19
%AoAhse0 is the tail stall angle when the elevator is not employed 
dAOAhE=4.2; %??????????????/
AOAhs=(AoAhse0-dAOAhE);
fprintf(" The tail stall angle of attack during take-off rotation (αhs)  =%8.3f deg  \n", AOAhs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Curve 

% % Most aft cg
xMgCgAf= Xcg;
xActCgAf=lm+Xcg;
xMgAc = Xacwf; % m from main landing gear
h_to_ho =(xMgAc-xMgCgAf)/mean_chord; % m
l_h1 = xActCgAf; %m
VH1 = (l_h1*Sh)/(Sw*mean_chord);
CmdE1 = -Clah*eh*VH1*te;
Cma1 = Clawf*h_to_ho-Clah*eh*Sh*(l_h1/mean_chord)*(1-de_da)/Sw;
% 
% Most forward cg
% xMgCgF= 0.5;
% xActCgF=1.3;
% xMgAc = 0.4; % m from main landing gear
% h_to_ho = (xMgAc-xMgCgF)/mean_chord; % m
% l_h2 = xActCgF; % m
% VH2 = (l_h2*Sh)/(Sw*mean_chord);
% CmdE2 = -Clah*eh*VH2*te;
% Cma2 = Clawf*h_to_ho-Clah*eh*Sh*(l_h2/mean_chord)*(1-de_da)/Sw;

i =1;


for U1=V_stall:0.1:Vmax
Q=0.5*ruo*U1^2;
CL1= (W)/(Q*Sw);
f1=((Tmax*0)/(Q*Sw*mean_chord))+Cm0; %##### Cm0
dE1(i)=-((f1*Clawf)+(CL1-Cl0)*Cma1)/(Clawf*CmdE1-Cma1*CldE); %##### Cl0
%dE2(i)=-((f1*Clawf)+(CL1-Cl0)*Cma2)/(Clawf*CmdE2-Cma2*CldE);

V(i)=U1;
i=i+1;
end
plot(V,dE1*57.3,'*')
grid
xlabel ('Speed (m/s)')
ylabel (' delta_E (deg)')
legend('Most aft cg')


%%%%%%%%%%%%5
N = 10; % (number of segments-1)
ARh = Sh/Bh^2; % Aspect ratio
%lambda = lambda; % Taper ratio #############
alpha_twist = -0.00001; % Twist angle (deg)
a_2d = 6.3; % lift curve slope (1/rad) ##############
a_0_e = 20; % def up zero-lift angle of attack (deg)
a_0_ed = -20; % def down zero-lift angle of attack (deg)
BE_BH=1; %ele-to-tail span ratio
Croot = (1.5*(1+lambda)*elevator_chord)/(1+lambda+lambda^2); % root chord
theta = pi/(2*N):pi/(2*N):pi/2;
alpha=ih+alpha_twist:-alpha_twist/(N-1):ih;
% segment’s angle of attack
for i=1:N
if (i/N)>(1-BE_BH)
alpha_0(i)=a_0_ed; %def down zero lift AOA
else
alpha_0(i)=a_0_e; %def up zero lift AOA
end
end
z = (elevator_span/2)*cos(theta);
c = Croot * (1 - (1-lambda)*cos(theta)); % MAC at each segment
mu = c * a_2d / (4 * elevator_span);
LHS = mu.*(alpha-alpha_0)/57.3; % Left Hand Side
% Solving N equations to find coefficients A(i):
for i=1:N
for j=1:N
B(i,j) = sin((2*j-1) * theta(i)) * (1 + (mu(i) *(2*j-1)) / sin(theta(i)));
end
end

A=B\transpose(LHS);
for i = 1:N
sum1(i) = 0;
sum2(i) = 0;
for j = 1 : N
sum1(i) = sum1(i) + (2*j-1) * A(j)*sin((2*j-1)*theta(i));
sum2(i) = sum2(i) + A(j)*sin((2*j-1)*theta(i));
end
end
CL_TO = pi * ARh * A(1)

end