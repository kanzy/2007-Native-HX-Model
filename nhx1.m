%% 2007-12-05 H->D NHX simulation & fitting
%% N(cl)<=>I(cl or op)<=>U(op) 

%% y1:N y2:I y3:U y4:H
%% k1:N->I k2:I->N k3:I->U k4:U->I

clear

global k1; 
global k2; 
global k3; 
global k4; 
global kch; 

hxTime=720; %unit: hr
temp=25     %unit: 'C
temp=temp+273.15    %unit: K
R=8.314;

%%initial values when [Deu]=0
k1_0=2e2;   
k2_0=1e5;   
k3_0=1e4;
k4_0=1e5;
kch=1e-2;

m_NI=1.25   %unit: kcal/mol/M
m_IU=1.20

for i=1:6
    Denat(i)=-0.3+0.3*i;    %[Denaturant]: [0 0.3 0.6 0.9 1.2 1.5]
    m1=1; m2=m1/exp(m_NI*4200*Denat(i)/(R*temp));
    m3=1; m4=m3/exp(m_IU*4200*Denat(i)/(R*temp));
    k1=m1*k1_0; k2=m2*k2_0; k3=m3*k3_0; k4=m4*k4_0;

%%calculate K & deltaG:
K_NI(i)=k1/k2 
deltaG_NI(i)=-R*temp*log(K_NI(i))/4200 %unit: kcal/mol
K_IU(i)=k3/k4
deltaG_IU(i)=-R*temp*log(K_IU(i))/4200 %unit: kcal/mol
deltaG_NU(i)=deltaG_NI(i) + deltaG_IU(i);


subplot(1,2,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Case1: I(cl) //the residue locate in more stable (closed in I) foldon

%%simulation using 'nhx101':
[t,y] = ode15s('nhx101',[0:300:hxTime*3600],[k2*k4/(k2*k4+k1*k4+k1*k3) k1*k4/(k2*k4+k1*k4+k1*k3) k1*k3/(k2*k4+k1*k4+k1*k3) 1]); 
% loglog(t/3600, y(:,1:3))  %plot N,I,U partition
% hold on
semilogx(t/3600,y(:,4),'b','LineWidth',3)
xlabel('Ex Time (hr)')
ylabel('%H')
hold on

%%single-exponential fitting:
sizer=size(y); sizer=sizer(1);
iniK=log(y(1,4)/y(sizer,4))/(t(sizer)-t(1));      
iniPara=[y(1,4); iniK];  
options = optimset('TolX', 1e-9, 'TolFun', 1e-15);
[fitPara,r1,r2,exitFlag,output]=lsqnonlin(@hxfit2, iniPara, [0;0],[], options, t, y(:,4));
fitA=fitPara(1);
fitk=fitPara(2);
k_ex_101(i)=fitk
%%calculate deltaG_ex:
deltaG_ex_101(i)=-R*temp*log(fitk/kch)/4200 %unit: kcal/mol
%%plot fitting curve:
semilogx(t/3600,fitA*exp(-fitk*t),'k','LineWidth',1)
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Case2: I(op) //the residue locate in the less stable (opened in I) foldon

%%simulation using 'nhx102':
[t,y] = ode15s('nhx102',[0:300:hxTime*3600],[k2*k4/(k2*k4+k1*k4+k1*k3) k1*k4/(k2*k4+k1*k4+k1*k3) k1*k3/(k2*k4+k1*k4+k1*k3) 1]); 
semilogx(t/3600,y(:,4),'g','LineWidth',3)
xlabel('Ex Time (hr)')
ylabel('%H')
hold on

%%single-exponential fitting:
sizer=size(y); sizer=sizer(1);
iniK=log(y(1,4)/y(sizer,4))/(t(sizer)-t(1));      
iniPara=[y(1,4); iniK];  
options = optimset('TolX', 1e-25, 'TolFun', 1e-25);
[fitPara,r1,r2,exitFlag,output]=lsqnonlin(@hxfit2, iniPara, [0;0],[], options, t, y(:,4));
fitA=fitPara(1);
fitk=fitPara(2);
k_ex_102(i)=fitk
%%calculate deltaG_ex:
deltaG_ex_102(i)=-R*temp*log(fitk/kch)/4200 %unit: kcal/mol
%%plot fitting curve:
semilogx(t/3600,fitA*exp(-fitk*t),'k','LineWidth',1)
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Case3: N(op) //the residue locate in the unprotected region

%%simulation using 'nhx103':
[t,y] = ode15s('nhx103',[0:300:hxTime*3600],[k2*k4/(k2*k4+k1*k4+k1*k3) k1*k4/(k2*k4+k1*k4+k1*k3) k1*k3/(k2*k4+k1*k4+k1*k3) 1]); 
semilogx(t/3600,y(:,4),'c','LineWidth',3)
xlabel('Ex Time (hr)')
ylabel('%H')
hold on

%%single-exponential fitting:
sizer=size(y); sizer=sizer(1);
iniK=log(y(1,4)/y(sizer,4))/(t(sizer)-t(1));      
iniPara=[y(1,4); iniK];  
options = optimset('TolX', 1e-25, 'TolFun', 1e-25);
[fitPara,r1,r2,exitFlag,output]=lsqnonlin(@hxfit2, iniPara, [0;0],[], options, t, y(:,4));
fitA=fitPara(1);
fitk=fitPara(2);
k_ex_103(i)=fitk
%%calculate deltaG_ex:
deltaG_ex_103(i)=-R*temp*log(fitk/kch)/4200 %unit: kcal/mol
%%plot fitting curve:
semilogx(t/3600,fitA*exp(-fitk*t),'k','LineWidth',1)
axis([0 1000 -0.05 1.05])
hold on

end

subplot(1,2,2)

plot(Denat, real(deltaG_ex_101),'ob')
hold on
plot(Denat, real(deltaG_ex_102),'og')
hold on
plot(Denat, real(deltaG_ex_103),'oc')
hold on
plot(Denat, deltaG_NI,'r')
hold on
plot(Denat, deltaG_NU,'b')
xlabel('[Denat] (M)')
ylabel('deltaG ex')




