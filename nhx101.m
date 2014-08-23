%% 2007-12-05 N(cl)<=>I(cl)<=>U(op) H->D NHX equations

function dy=nhx101(t,y)
dy = zeros(4,1);
%% y1:N y2:I y3:U y4:H

global k1; 
global k2; 
global k3; 
global k4; 
global kch; 

dy(1) = k2*y(2) - k1*y(1);
dy(2) = k1*y(1) + k4*y(3) - (k3+k2)*y(2);
dy(3) = k3*y(2) - (k4+kch)*y(3);
dy(4) = -kch*y(3);