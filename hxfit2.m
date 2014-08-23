%%hxfit2.m 2007-12-05

function f = hxfit2(x, time, data)

A=x(1);
K=x(2);

funxdata = A*exp(-K*time);

% sizer=size(data); sizer=sizer(1);
% for i=1:sizer
%     if data(i)<0
%         data(i)=1e-8;
%     end
% end
     
   f = data - funxdata;