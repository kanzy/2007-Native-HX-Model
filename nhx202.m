function dy=nhx202(t,y)
dy = zeros(5,1);
%% y1:N y2:I y3:U y4:H y5:M

global k1; 
global k2; 
global k3; 
global k4; 
global k5;
global k6;
global kch; 

dy(1) = k2*y(2) - k1*y(1);
dy(2) = k1*y(1) + k4*y(3) - (k3+k2+k5)*y(2) + k6*y(5);
dy(3) = k3*y(2) - (k4+kch)*y(3);
dy(4) = -kch*(y(3)+y(5));

dy(5) = k5*y(2) - (k6+kch)*y(5);