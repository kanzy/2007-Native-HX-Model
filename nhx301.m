function dy=nhx301(t,y)
dy = zeros(6,1);
%% y1:N y2:I y3:U y4:H y5:N' y6:I'

global k1; 
global k2; 
global k3; 
global k4; 
global k7;
global k8;
global k9;
global k10;
global kch; 

dy(1) = k2*y(2) - (k1+k7)*y(1) + k8*y(5);
dy(2) = k1*y(1) + k4*y(3) - (k3+k2+k9)*y(2) + k10*y(6);
dy(3) = k3*y(2) - (k4+kch)*y(3);
dy(4) = -kch*(y(3)+y(5)+y(6));

dy(5) = k7*y(1) - (k8+kch)*y(5);
dy(6) = k9*y(2) - (k10+kch)*y(6);