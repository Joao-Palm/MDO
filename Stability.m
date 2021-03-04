function [Max_stab,ME_min,Cl_trim,alpha_limite]=Stability(wing,stabilizer,GMP,downwash,alpha_0w);
         
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Facilitador de variáveis
z_h=stabilizer.geo.ref_point(3);
z_w=wing.geo.ref_point(3);
x_h=stabilizer.geo.ref_point(1);
x_w=wing.geo.ref_point(1);
b_w=wing.ref.b_ref;
Cr_w=wing.geo.c;
Ct_w=Cr_w*wing.geo.T(1,2);
S_w=wing.ref.S_ref;
S_h=stabilizer.ref.S_ref;
Sweep_w=wing.geo.SW(1,1);
alpha=wing.polar{1,1}(:,1);
v=wing.state.AS;
rho=wing.state.rho;
q=0.5*rho*v^2;
ip=GMP.ip*pi/180;
AR_w=(b_w^2)/S_w;
auxDist = abs(GMP.posx - wing.geo.ref_point(1))/wing.ref.C_mac;
alpha_min_voo=-3;

CA_stabilizer=stabilizer.geo.startx+stabilizer.geo.c/4;
CA_wing=wing.geo.startx+Cr_w/4;
%alpha_0w=interp1(wing.polar{1,1}(:,2),wing.polar{1,1}(:,1),0,'pchip');

%% Calculo do upwash   %top
upwash_a = upwashDerivative(AR_w, auxDist);
upwash = upwash_a * (alpha - (alpha_0w)); %Angulos estão em radianos

%% GMP %top
dxp=(GMP.posx-wing.geo.CG(1))/wing.ref.C_mac; %top
dzp=(GMP.posz-wing.geo.CG(3))/wing.ref.C_mac; %top

Tp=GMP.T_stat+GMP.d_t1*v; %top
Kt=Tp/(rho*(GMP.diam_p*v)^2); %top
if Kt>2.5
    Kt=2.5;
end
ftp=0.0442*Kt^3-0.2609*Kt^2+0.7833*Kt+1.0253; %top
Cn_blade=0.00094*GMP.jp^3-0.01191*GMP.jp^2+0.05567*GMP.jp; %top
Sp=(pi*GMP.diam_p^2)/4;
fpa=q*GMP.n_hel*Sp*Cn_blade*ftp; %top'
fp=fpa*(ip + alpha+upwash);%top
CMCG_gmp=GMP.n_motor*(dxp*(-Tp*sin(ip)-fp*cos(ip))+dzp*(-Tp*cos(ip)+fp*sin(ip)))/(S_w*q); %top

%% Stabilizer top
alpha_ef=alpha-downwash;

dxh=(CA_stabilizer-wing.geo.CG(1))/wing.ref.C_mac;
dzh=(stabilizer.geo.startz-wing.geo.CG(3))/wing.ref.C_mac;
eta_h=max(1-(2.42*(wing.polar{1,1}(:,4).^0.5)/((CA_stabilizer-CA_wing-0.75*wing.ref.C_mac)/Cr_w)+0.3),0.85);
for i=1:stabilizer.n_deflex 
    for j=1:length(stabilizer.polar{1,i})
        cm1_h(j,i)=dxh*(-stabilizer.polar{1,i}(j,2)*cos(alpha_ef(j,1))-(stabilizer.polar{1,i}(j,3)+stabilizer.polar{1,i}(j,4))*sin(alpha_ef(j,1))); 
        cm2_h(j,i)=dzh*((stabilizer.polar{1,i}(j,3)+stabilizer.polar{1,i}(j,4))*cos(alpha_ef(j,1))-stabilizer.polar{1,i}(j,2)*sin(alpha_ef(j,1)));
        CMCG_h(j,i)=(cm1_h(j,i)+cm2_h(j,i)+stabilizer.polar{1,i}(j,5)*stabilizer.ref.C_mac/wing.ref.C_mac).*eta_h(j)*S_h/S_w;
    end
end
%% Wing top

dxw=(CA_wing-wing.geo.CG(1))/wing.ref.C_mac;
dzw=(wing.geo.startz-wing.geo.CG(3))/wing.ref.C_mac;

cm1_w=dxw*(-wing.polar{1,1}(:,2).*cos(alpha)-(wing.polar{1,1}(:,3)+wing.polar{1,1}(:,4)).*sin(alpha));
cm2_w=dzw*(-wing.polar{1,1}(:,2).*sin(alpha)+(wing.polar{1,1}(:,3)+wing.polar{1,1}(:,4)).*cos(alpha));
CMCG_w=cm1_w+cm2_w+wing.polar{1,1}(:,5);

%% Aircraft 

%Polar do avião - alpha x deflexão
 for j=1:stabilizer.n_deflex
        alpha_limite(j)=length(stabilizer.polar{1,j});
 end    
alpha_limite=min(alpha_limite);

for i=1:stabilizer.n_deflex
    for j=1:alpha_limite
        Cm_aircraft(j,i)=CMCG_w(j,1)+CMCG_h(j,i)+CMCG_gmp(j,1);
        Cl_aircraft(j,i)=wing.polar{1,1}(j,2)+stabilizer.polar{1,i}(j,2).*eta_h(j,1)*S_h/S_w;%+Tp*sin(alpha(j,1)+ip)/(0.5*rho*v^2*S_w);
        Cd_aircraft(j,i)=(wing.polar{1,1}(j,3)+wing.polar{1,1}(j,4))+(stabilizer.polar{1,i}(j,3)+stabilizer.polar{1,i}(j,4)).*eta_h(j,1)*S_h/S_w+Tp*cos(alpha(j,1)+ip)/(0.5*rho*v^2*S_w);
    end
end

for j=1:stabilizer.n_deflex
    for i=2:alpha_limite-1
            ME_matrix(i-1,j)=-100*(Cm_aircraft(i+1,j)-Cm_aircraft(i-1,j))/(Cl_aircraft(i+1,j)-Cl_aircraft(i-1,j));
    end
end


% Polar trimada - alpha x deflexão
for i=1:alpha_limite
    j=1;
    Trim(i)=0;
    while Trim(i)==0 && j<=stabilizer.n_deflex-1
        if Cm_aircraft(i,j)*Cm_aircraft(i,j+1)<=0
            Deflex(i,1)=interp1([Cm_aircraft(i,j) Cm_aircraft(i,j+1)],[j j+1],0);
            Clh_trim(i,1)=interp1([j j+1],[stabilizer.polar{1,j}(i,2) stabilizer.polar{1,j+1}(i,2)],(Deflex(i))); %CL do profundor p/ CM=0 
            Cdh_trim(i,1)=interp1([j j+1],[stabilizer.polar{1,j}(i,3) stabilizer.polar{1,j+1}(i,3)],(Deflex(i)))+interp1([j j+1],[stabilizer.polar{1,j}(i,4) stabilizer.polar{1,j+1}(i,4)],(Deflex(i))); %CD do profundor p/ CM=0
            Cl_trim(i,1)=wing.polar{1,1}(i,2)+(Clh_trim(i,1).*cos(downwash(i,1))-Cdh_trim(i,1).*sin(downwash(i,1))).*eta_h(i)*S_h/S_w;
            Cd_trim(i,1)=(wing.polar{1,1}(i,3)+wing.polar{1,1}(i,4))+Cdh_trim(i,1).*eta_h(i)*S_h/S_w+Tp*cos(alpha(i,1)+ip)/(0.5*rho*v^2*S_w);
            Cl_trim(i,2)=alpha(i,1)*180/pi;
            Stabilizer_use(i,1)=abs(100*Clh_trim(i,1)./max(abs(stabilizer.polar{1,j}(:,2))));
            if i<=length(ME_matrix)
                ME(i,1)=interp1([j j+1],[ME_matrix(i,j) ME_matrix(i,j+1)],Deflex(i));
            end
            Trim(i)=1;
        else
            Deflex(i,1)=NaN;
%             Clh_trim(i,1)=NaN;
%             Cdh_trim(i,1)=NaN;
%             Stabilizer_use(i,1)=NaN;
%             Cl_trim(i,1)=NaN;
%             ME(i,1)=NaN;
        end
        j=j+1;
    end
    
end

%Margem Estática [-]

ME_min=interp1(Cl_trim(:,1),ME(:,1),0.5,'linear','extrap');

pos=find(wing.polar{1,1}(:,1) == 0);
Max_stab=max(Stabilizer_use((pos+alpha_min_voo):end));
