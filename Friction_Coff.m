function [Cf,Swett] = Friction_Coff(y,wing_span,Twist_angel,Alpha_induce_0,WSSL,density,cord,v_cruise,air_viscosity,i_w,Alpha_FOP,diameter,high_wing,RC_w,sign)
pi_Twist = (2.*y./wing_span).*Twist_angel;                                  % span wist twist relative to root
Alpha_induce = flip(Alpha_induce_0);                                        % induce angle of attack along semi span
WSSL = flip(WSSL);                                                          % semi span induce angle of attack locations
[row col] = size(y);
if col>2
    for i = 2:(col-1)
        larger_Loc = find(WSSL>y(i));
        larger_Loc = min(larger_Loc);
        smaller_Loc = find(WSSL<y(i));
        smaller_Loc = max(smaller_Loc);
        Larger_ai = Alpha_induce(larger_Loc);
        smaller_ai = Alpha_induce(smaller_Loc);
        ai(i-1) = (((y(i)-WSSL(smaller_Loc))*(Larger_ai-smaller_ai))/...
            (WSSL(larger_Loc)-WSSL(smaller_Loc))) + smaller_ai;
    end
    Induced_AOA = [Alpha_induce(1) ai Alpha_induce(end)];
else
    Induced_AOA = [Alpha_induce(1) Alpha_induce(end)];
end
%%
Re =  (density).*cord.*v_cruise ./ air_viscosity;                           % Rynold number @ Root Airfoil
Eff_A = i_w+pi_Twist-Alpha_FOP-Induced_AOA;                                 % effective angle of attack @ root airfoil


for i = 1:col
    if sign == 1
        j = num2str(1);
    elseif sign == 2
        j = num2str(2);
    elseif sign == 3
        j = num2str(3);
    end

        saveFlnAF_data = num2str(i);
    fidAF = fopen([saveFlnAF_data j '.txt']);                                     % Open file for reading
    dataBuffer = textscan(fidAF,'%f %f %f %f %f %f %f ' ,'HeaderLines',12,...   % open data file
        'CollectOutput',1,...
        'Delimiter','');
    fclose(fidAF);                                                              % Close the file
    
    alpha  = dataBuffer{1,1}(:,1);                                              % angle of attack
    XT_Top = dataBuffer{1,1}(:,6);                                              % transition point on upper surface
    XT_Bottom = dataBuffer{1,1}(:,7);                                           % transition point on lower surface

    larger_Loc = find(alpha>Eff_A(i));
    
    if Eff_A(i)>=0
        larger_Loc = min(larger_Loc);
    else
        larger_Loc = max(larger_Loc);
    end
    
    smaller_Loc = find(alpha<Eff_A(i));
    
    if Eff_A(i)>=0
        smaller_Loc = max(smaller_Loc);
    else
        smaller_Loc = min(smaller_Loc);
    end
    
    Larger_XT = XT_Top(larger_Loc);
    smaller_XT = XT_Top(smaller_Loc);
    XT_Top = (((Eff_A(i)-alpha(smaller_Loc))*(Larger_XT-smaller_XT))/...
        (alpha(larger_Loc)-alpha(smaller_Loc))) + smaller_XT;
    
    Larger_XB = XT_Bottom(larger_Loc);
    smaller_XB = XT_Bottom(smaller_Loc);
    XT_Bottom = (((Eff_A(i)-alpha(smaller_Loc))*(Larger_XB-smaller_XB))/...
        (alpha(larger_Loc)-alpha(smaller_Loc))) + smaller_XB;
    
    
    
    X0_CT = 36.9*(XT_Top)^0.625 *(1/Re(i))^0.375;                               % Fictitious Turbulent BL on Root Airfoil  Upper Surface
    X0_CB = 36.9*(XT_Bottom)^0.625 *(1/Re(i))^0.375;                            % Fictitious Turbulent BL on Root Airfoil  Lower Surface
    Cf_Uc = (0.074/Re(i)^0.2)*(1-((XT_Top)-X0_CT))^0.8;                         % Skin Friction for Root Airfoil  Upper Surface
    Cf_Bc = (0.074/Re(i)^0.2)*(1-((XT_Bottom)-X0_CB))^0.8;                      % Skin Friction for Root Airfoil  Lower Surface
    Cf_c(i) = 0.5*(Cf_Uc+Cf_Bc);                                                % Average Skin Friction for Root Airfoil
end


[ro,co] = size(cord);
for i = 1:(co-1)                                                            % from 1:number of segmens 
    CFi(i) =0.5*(Cf_c(i)+Cf_c(i+1));
    Sref(i) = 0.5*(cord(i)+cord(i+1))*(y(i+1)-y(i))*2;
    Swett(i) = 2*Sref(i);
end

if high_wing ==1
    Swett(1) = Swett(1)-(RC_w*diameter);
else
    Swett(1) = Swett(1)-(RC_w*diameter*2);
end
Swett = 1.1*Swett;
Cf_sw = CFi.*Swett;
Swett = sum(Swett);
Cf = sum(Cf_sw)/sum(Swett);
end