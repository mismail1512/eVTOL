
function [CL_wing,e,Alpha_induce,z] = Linear_Lifting_Line_Method(AR_w,RC_w,wing_span,lambda,i_w,Twist_angel,a_2d,alpha_0,yy,Alpha_FOP,density,velocity)
N = 30; % (number of segments - 1)
theta = pi/(2*N):pi/(2*N):pi/2;
alpha = i_w-Alpha_FOP+Twist_angel:-Twist_angel/(N-1):i_w-Alpha_FOP;
% segment’s angle of attack
z = (wing_span/2)*cos(theta);
c = RC_w * (1 - (1-lambda)*cos(theta)); % Mean AerodynamicsChord at each segment (m)
mu = c * a_2d / (4 * wing_span);
LHS = mu .* (alpha-alpha_0)/57.3; % Left Hand Side
% Solving N equations to find coefficients A(i):
for i=1:N
    for j=1:N
        B(i,j) = sin((2*j-1) * theta(i)) * (1 + (mu(i) * (2*j-1)) /sin(theta(i)));
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
CL = 4*wing_span*sum2 ./ c;

% span efficiency factor
se = 0;
for i = 3:2:N
    
    se = se + i*(A(i)/A(1));
end
e = 1/(1+se);

% induce angle of attack 
for i = 1:N
    induce_angle = 0;
    for j = 1:N
        induce_angle = induce_angle+ j*A(j)*sin(j*theta(i))./sin(theta(i));
    end
    Alpha_induce(i) = sum(induce_angle).*57.3;
end
w = -Alpha_induce.*velocity;

%% Plotting 
CL1 =[0];
y_s = [wing_span/2];
for ii = 1:N
    CL1 = [CL1,CL(ii)];
    y_s = [y_s,z(ii)];
end
c = [0 c];
L = 0.5*density*velocity^2.*c.*CL1;
f = -y_s;
f = flip(f);
%y_s = [f y_s];
l = flip(L);
L = [l L];
o = flip(CL1);
%CL1 = [o CL1];
if(yy ==1)
    figure
    hold all
    plot(y_s,CL1,'-o')
    grid
    title('Lift distribution')
    xlabel('span location (m)')
    ylabel ('Lift Coefficient')
end
CL_wing = pi * AR_w * A(1);
end