%% y1:N y2:I y3:U y4:H y5:N' y6:I' y7:M
%% k1:N->I k2:I->N k3:I->U k4:U->I k5:I->M k6:M->I k7:N->N' k8:N'->N
%% k9:I->I' k10:I'->I

clear

global k1; 
global k2; 
global k3; 
global k4; 
global k5;
global k6

global k7;
global k8;
global k9;
global k10;
global kch; 

%%residues=[# k7   k8   k9   k10 case]
residues=[1  0e2  1e6  0e2  1e6  1  
          2  1e1  1e6  1e1  1e6  1
          3  1e3  1e6  1e3  1e6  1
          4  1e5  1e6  1e5  1e6  1
          5  1e2  0e6  1e2  0e6  1
          6  0e2  1e6  0e2  1e6  2  
          7  1e1  1e6  1e1  1e6  2
          8  1e3  1e6  1e3  1e6  2
          9  2e5  1e6  2e5  1e6  2
          10 1e2  0e6  1e2  0e6  2
          11 0e2  1e6  0e2  0e6  3  
          12 1e1  1e6  0e1  0e6  3
          13 1e3  1e6  0e3  0e6  3
          14 2e5  1e6  0e5  0e6  3
          15 1e2  0e6  0e2  0e6  3
          16 0e2  1e6  0e2  0e6  4  
          17 1e1  1e6  0e1  0e6  4
          18 1e3  1e6  0e3  0e6  4
          19 3e5  1e6  0e5  0e6  4
          20 1e2  0e6  0e2  0e6  4
      ];

Denat=[0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4 5];    %[Denaturant] unit:M 

hxTime=1200; %unit: hr
temp=25;     %unit: 'C
temp=temp+273.15;    %unit: K
R=8.314;

kch=1e-2;

k1_0=1e2;   %N->I unfolding
k2_0=1e6;   %I->N folding
k3_0=1e3;   %I->U unfolding
k4_0=1e6;   %U->I folding

k6_0=1e3;   %M->I unfolding
k5_0=1e4;   %I->M folding

m_NI=1.25;   %unit: kcal/mol/M
m_IU=1.20;
m_IM=1.2;

subplot(1,2,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%analysis loop for different [Denaturant]:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sizer=size(Denat); DenatPoints=sizer(2);
for i=1:DenatPoints
    
    m1=1; m2=m1/exp(m_NI*4200*Denat(i)/(R*temp));
    m3=1; m4=m3/exp(m_IU*4200*Denat(i)/(R*temp));
    m6=1; m5=m6/exp(m_IM*4200*Denat(i)/(R*temp));
    k1=m1*k1_0; k2=m2*k2_0; k3=m3*k3_0; k4=m4*k4_0; k5=m5*k5_0; k6=m6*k6_0;
    
    N0=k2*k4*k6/(k1*k4*k6+k2*k4*k6+k1*k3*k6+k1*k4*k5)
    I0=k1*k4*k6/(k1*k4*k6+k2*k4*k6+k1*k3*k6+k1*k4*k5)
    M0=k1*k4*k5/(k1*k4*k6+k2*k4*k6+k1*k3*k6+k1*k4*k5)
    U0=k1*k3*k6/(k1*k4*k6+k2*k4*k6+k1*k3*k6+k1*k4*k5)
        
    %%calculate K & deltaG:
    K_NI(i)=k1/k2; 
    deltaG_NI(i)=-R*temp*log(K_NI(i))/4200; %unit: kcal
    K_IU(i)=k3/k4;
    deltaG_IU(i)=-R*temp*log(K_IU(i))/4200; %unit: kcal
    deltaG_NU(i)=deltaG_NI(i) + deltaG_IU(i);
    K_IM(i)=k6/k5;
    deltaG_IM(i)=-R*temp*log(K_IM(i))/4200; %unit: kcal
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%analysis loop for different residues: 
    sizer=size(residues); residueNum=sizer(1);
    for j=1:residueNum
        
        k7=residues(j,2);
        k8=residues(j,3);
        k9=residues(j,4);
        k10=residues(j,5);
        
        residue_case=residues(j,6);
        
        %%simulation:
        switch residue_case
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%Case1: I(cl)&M(cl) //the residue locate in
            case 1
                [t,y] = ode15s('nhx401',[0:600:hxTime*3600],[N0 I0 U0 1 0 0 M0]); 
                semilogx(t/3600,y(:,4),'b','LineWidth',2)
                hold on
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%Case2: I(cl)&M(op) //the residue locate in
            case 2
                [t,y] = ode15s('nhx402',[0:600:hxTime*3600],[N0 I0 U0 1 0 0 M0]); 
                semilogx(t/3600,y(:,4),'g','LineWidth',2)
                hold on
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%Case3: I(op)&M(cl) //the residue locate in
            case 3
                [t,y] = ode15s('nhx403',[0:600:hxTime*3600],[N0 I0 U0 1 0 0 M0]); 
                semilogx(t/3600,y(:,4),'c','LineWidth',2)
                hold on
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%Case4: I(op)&M(op) //the residue locate in
            case 4
                [t,y] = ode15s('nhx404',[0:600:hxTime*3600],[N0 I0 U0 1 0 0 M0]); 
                semilogx(t/3600,y(:,4),'y','LineWidth',2)
                hold on
        end
        
        %%single-exponential fitting:
        sizer=size(y); dataSize=sizer(1);
        iniK=log(y(1,4)/y(dataSize,4))/(t(dataSize)-t(1));      
        iniPara=[y(1,4); iniK];  
        options = optimset('TolX', 1e-9, 'TolFun', 1e-15);
        [fitPara,r1,r2,exitFlag,output]=lsqnonlin(@hxfit2, iniPara, [0;0],[], options, t, y(:,4));
        fitA=fitPara(1);
        fitk=fitPara(2);
        k_ex(i,j)=real(fitk);
        %%calculate deltaG_ex:
        deltaG_ex(i,j)=-R*temp*log(fitk/kch)/4200; %unit: kcal/mol
        %%plot fitting curve:
        semilogx(t/3600,fitA*exp(-fitk*t),'k','LineWidth',1)
        hold on
    i
    j    
    end
    
end
axis([0 1500 -0.05 1.05])
xlabel('Ex Time (hr)')
ylabel('%H')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2)

plot(Denat, deltaG_NI,'r','LineWidth',2)
hold on
plot(Denat, deltaG_NU,'k','LineWidth',2)
hold on
plot(Denat, deltaG_IM,'y','LineWidth',2)

for j=1:residueNum
    switch residues(j,6)
        case 1
            plot(Denat, real(deltaG_ex(:,j)),'--ob','LineWidth',1)
        case 2
            plot(Denat, real(deltaG_ex(:,j)),'--^b','LineWidth',1)
        case 3
            plot(Denat, real(deltaG_ex(:,j)),'--om','LineWidth',1)
        case 4
            plot(Denat, real(deltaG_ex(:,j)),'--^m','LineWidth',1)
    end
    hold on
end

xlabel('[Denat] (M)')
ylabel('deltaG ex (kcal/mol)')