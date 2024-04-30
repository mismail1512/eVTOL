
function [t2,b_A,C_A, Aa,yi,y0] = Aileronn(S_w,lambda,S_Ht,S_Vt,v_stall,cl_alpha_3d_wf,...
    I_xx,wing_span,Cr_w,Phai_des,CWing)
fprintf("***********Aileron **************\n");

plott=1;
Taw=0.4;
bai_b=0.65;
bao_b=0.95;
Aileron_defle = 25;
Ca_C = 0.2;
yi = bai_b*wing_span/2;                                  % Inboard location of aileron span
y0 = bao_b*wing_span/2;                                  % Outboard Location of Aileron span 

CL_A = (2*cl_alpha_3d_wf*Taw*Cr_w/(wing_span*S_w))*((y0^2 /2 +(2/3)*y0^3 *((lambda-1)/wing_span))...
-(yi^2 /2 +(2/3)*yi^3 *((lambda-1)/wing_span)));           % Aileron rolling Moment coff derivative with respect to aileron deflection

CL = CL_A*(Aileron_defle/57.3);                     % Rolling Moment Coff
V_app = 1.3*v_stall;                              % approach speed
L_A = 0.5*1.1117*(V_app^2) *S_w*CL*wing_span;         % Rolling Moment due to Aileron Max def
fprintf("The aircraft rolling moment (LA) when the aileron is deflected with the maximum deflection is: (=%8.3f ) \n", L_A);


Y_D =0.6*wing_span/2;    %###################                               % Drag Moment Arm
CD_R = 0.8;  %?????????
P_ss = sqrt(2*L_A/(1.1117*(S_w+S_Ht+S_Vt)*CD_R*Y_D^3));      % Steady State Roll Rate(Rad/sec)
Phai1 = log(P_ss^2)*I_xx/(1.1117*(S_w+S_Ht+S_Vt)*CD_R*Y_D^3);     % Bank Angle @ which Aircraft achieve steady state Roll Rate(rad)
fprintf("The bank angle (1) at which the aircraft achieves the steady-state roll rate is: (=%8.3f deg) \n", Phai1);


P_dot = P_ss^2 /(2*Phai1);                 % aircraft angular acceleration before steadt state roll
fprintf("the aircraft rate of roll rate (•P) that is produced by the aileron rolling moment until the aircraft reaches the steady-state roll rate (P ss) is: (=%8.3f deg^2) \n", P_dot);
n=rem(Phai1*57.3, 360 );
if Phai1*57.3 > Phai_des
    t2 = sqrt(2*(Phai_des/57.3)/P_dot);                 % time taken to achieve desired bank angle 
    fprintf("he time it takes the aircraft to achieve the bank angle of (=%8.3f deg) deg is (=%8.3f ) s \n", Phai_des,t2);

else
    t1 = sqrt(2*(Phai1)/P_dot);
    dt21 = (Phai_des/57.3-Phai1)/P_ss;
    t2 = t1+dt21;
    t= sqrt(2*(Phai_des)/P_dot);
    fprintf("he time it takes the aircraft to achieve the bank angle of (=%8.3f deg) deg is (=%8.3f ) s \n", Phai_des,t2);

end

if plott ==1
figure(4)
PI_des = linspace(0,40,27);
T2 = sqrt(2.*(PI_des./57.3)./P_dot);
plot(T2,PI_des,'*');
grid;
ylabel('Bank angle');
xlabel ('time');
end
b_A = y0-yi;
fprintf("The span of ailerons is: (=%8.3f)m \n", b_A);
C_A = Ca_C*CWing;
fprintf("The chord of ailerons is: (=%8.3f)m \n", C_A);
Aa=2*b_A*C_A;
fprintf("The overall planform area of both left and right ailerons is: (=%8.3f m) \n", Aa);

end


        