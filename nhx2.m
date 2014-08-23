%% y1:N y2:I y3:U y4:H y5:M
%% k1:N->I k2:I->N k3:I->U k4:U->I k5:I->M k6:M->I

clear

global k1; 
global k2; 
global k3; 
global k4; 

global k5;
global k6;
global kch; 

hxTime=720; %unit: hr
temp=25;     %unit: 'C
temp=temp+273.15;    %unit: K
R=8.314;

k1_0=2e2;   
k2_0=1e5;   
k3_0=1e4;
k4_0=1e5;

k5_0=1e0;
k6_0=3e1;
kch=1e-2;

m_NI=1.25;   %unit: kcal/mol/M
m_IU=1.20;
m_IM=1.2;

for i=1:6
    Denat(i)=-0.3+0.3*i;    %[Denaturant]: [0 0.3 0.6 0.9 1.2 1.5]
    m1=1; m2=m1/exp(m_NI*4200*Denat(i)/(R*temp));
    m3=1; m4=m3/exp(m_IU*4200*Denat(i)/(R*temp));
    m6=1; m5=m6/exp(m_IM*4200*Denat(i)/(R*temp));
    k1=m1*k1_0; k2=m2*k2_0; k3=m3*k3_0; k4=m4*k4_0; k5=m5*k5_0; k6=m6*k6_0;

%%calculate K & deltaG:
K_NI(i)=k1/k2;
deltaG_NI(i)=-R*temp*log(K_NI(i))/4200; %unit: kcal
K_IU(i)=k3/k4;
deltaG_IU(i)=-R*temp*log(K_IU(i))/4200; %unit: kcal
K_IM(i)=k6/k5;
deltaG_IM(i)=-R*temp*log(K_IM(i))/4200; %unit: kcal
deltaG_NU(i)=deltaG_NI(i) + deltaG_IU(i);

subplot(1,2,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Case1: I(cl)&M(cl) //the residue locate in

%%simulation using 'nhx201':
[t,y] = ode15s('nhx201',[0:300:hxTime*3600],[1 0 0 1 0]); 
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
k_ex_201(i)=fitk
%%calculate deltaG_ex:
deltaG_ex_201(i)=-R*temp*log(fitk/kch)/4200 %unit: kcal/mol
%%plot fitting curve:
semilogx(t/3600,fitA*exp(-fitk*t),'k','LineWidth',1)
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Case2: I(cl)&M(op) //the residue locate in

%%simulation using 'nhx202':
[t,y] = ode15s('nhx202',[0:300:hxTime*3600],[1 0 0 1 0]); 
semilogx(t/3600,y(:,4),'g','LineWidth',3)
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
k_ex_202(i)=fitk
%%calculate deltaG_ex:
deltaG_ex_202(i)=-R*temp*log(fitk/kch)/4200 %unit: kcal/mol
%%plot fitting curve:
semilogx(t/3600,fitA*exp(-fitk*t),'k','LineWidth',1)
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Case3: I(op)&M(cl) //the residue locate in

%%simulation using 'nhx203':
[t,y] = ode15s('nhx203',[0:300:hxTime*3600],[1 0 0 1 0]); 
semilogx(t/3600,y(:,4),'c','LineWidth',3)
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
k_ex_203(i)=fitk
%%calculate deltaG_ex:
deltaG_ex_203(i)=-R*temp*log(fitk/kch)/4200 %unit: kcal/mol
%%plot fitting curve:
semilogx(t/3600,fitA*exp(-fitk*t),'k','LineWidth',1)
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Case4: I(op)&M(op) //the residue locate in

%%simulation using 'nhx204':
[t,y] = ode15s('nhx204',[0:300:hxTime*3600],[1 0 0 1 0]); 
semilogx(t/3600,y(:,4),'m','LineWidth',3)
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
k_ex_204(i)=fitk
%%calculate deltaG_ex:
deltaG_ex_204(i)=-R*temp*log(fitk/kch)/4200 %unit: kcal/mol
%%plot fitting curve:
semilogx(t/3600,fitA*exp(-fitk*t),'k','LineWidth',1)
hold on
axis([0 1000 -0.05 1.05])

end

subplot(1,2,2)

plot(Denat, real(deltaG_ex_201),'ob')
hold on
plot(Denat, real(deltaG_ex_202),'og')
hold on
plot(Denat, real(deltaG_ex_203),'oc')
hold on
plot(Denat, real(deltaG_ex_204),'om')
hold on
plot(Denat, deltaG_NI,'r')
hold on
plot(Denat, deltaG_NU,'b')
xlabel('[Denat] (M)')
ylabel('deltaG ex')