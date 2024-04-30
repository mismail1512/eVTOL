function [D_P,L_P,S_p,D_Cyl,L_Cyl,S_Cyl,D_CO,L_CO,S_CO,D_fL,D_fs,L_f,S_f]...
    = Basic_Geom(R2,R1,length,RC_w,X_cg_dash,C_MGC,lt)
% parabloid 
D_P = 2*R2;                                                                 % parabloid base diameter in (m)
L_P = 0.28*length;                                                          % parabloid height (m)
S_p =((pi*D_P/(12*L_P^2))*((4*L_P^2 + D_P^2/4)^1.5-D_P^3/8)-pi*D_P^2/4);    % parabloid lateral surface area
% circular cylinder  
D_Cyl = 2*R2;                                                               % circular cylinder diameter (m)
L_Cyl = RC_w;                                                               % circular cylinder heigh (m)
S_Cyl = pi*D_Cyl*L_Cyl;                                                     % cylinder lateral surface area                                                       
% cone 
LL = (lt-(RC_w-(X_cg_dash*C_MGC)));
D_CO = 0.30;                                                                % cone base diameter(m)
L_CO = LL;                                                                  % cone height
S_CO = (pi*D_CO/2)*sqrt(L_CO^2 +D_CO^2/4);                                  % cone lateral surface area 
% Frustum 
D_fL = [2*R2 R2];                                                           % [ left frustum_larger radius , right frustum_larger diameter] 
D_fs = [R2 R1];                                                             % [ left frustum_smaller radius ,right frustum_smaller diameter] 
L_f = [(1/3)*LL (2/3)*LL];                   
S_f = (pi*(D_fL + D_fs)/2).*sqrt(L_f.^2 + (D_fL.^2 -D_fs.^2)./4);           % [left frustum area,right frustum area]