function [Gear_hight, z_nose_gear , z_main_gear] = designLandingGear( tail_arm )

% the gear configuration is choosen to be fixed tricycle .. 
% assuming the take off angle to be in average 10 degrees
% so .. 
takeoff_angle = 10 ;   % degrees
Gear_hight = tail_arm * sin(takeoff_angle * pi / 180);
% assume the tail arm to the fueslage length ratio to be 0.65 
fueslage_length = tail_arm / 0.65 ; 
% let the nose gear to be at 5 % of the fueselage length from the nose 
% and the Cg is assumed to be at 35 % of the fueselage length

z_nose_gear = 0.05 * fueslage_length;      % the distance between the cg and the nose gear. 

% the nose gear should carry aroung 15 % of the aircraft weight .. 
z_main_gear = ( ( 0.15 * 0.3 ) / 0.85 ) * fueslage_length + 0.35 * fueslage_length;  

fprintf("The Fueselage length  = %f m\n",fueslage_length);
fprintf("The nose Gear location at %f m from the nose \n", z_nose_gear );
fprintf("The main Gear location at %f m from the nose \n",z_main_gear);
fprintf("The Gear Hight is %f m \n \n",Gear_hight);

end