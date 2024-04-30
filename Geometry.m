function [span,C_avg,RC,TC,C_MGC] = Geometry(AR,S,tr)
span = (AR*S)^0.5;                                             % horizontal tail span 
C_avg = S/span;                                                % horizontal tail average cord
RC=2*C_avg/(1+tr);                                                 % horizontal tail root cord 
TC=RC*tr;                                                          % horizontal tail tip cord                                                         
C_MGC=(2/3)*RC*(tr^2+tr+1)/(tr+1);
end   