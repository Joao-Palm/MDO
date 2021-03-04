 clear
% clc
%% MDO - Alpha 2020 - Main Code
%Author: João Gonçalves


%Comentários importantes: -Se mudar o nº de painéis, modificar o taper ratio no downwash 
%% Inputs
cd('D:\AeroDesign\ALPHA 2020\MDO 2020')
alpha_inicial= -5;
alpha_incremento=1;
alpha_final=25;



xCG=linspace(0.005,0.045,60);
zCG=linspace(0.05,0.105,10);

% xCG=[0.009];
% zCG=[0.1];

constrains.ME=5; % [%]
constrains.Stab=80; %[%]

airfoil_lib(1)={'S21 pomba'};
airfoil_lib(2)={'maggie watson pomba'};
% airfoil_lib(3)={'Watson'};
% airfoil_lib(4)={'Zev'};
% airfoil_lib(5)={'Abapa51'};
% airfoil_lib(6)={'Alpha013'};
% airfoil_lib(7)={'SigFRA'};
airfoil_lib(3)={'SigLHR'};
airfoil_lib(4)={'SigPVG'};

%modeFrontier Wing

Tr1_w=0.95;
Tr2_w=0.92;
b_h=1.62;
b_w=2.19;
cr_h=0.28;
cr_w=0.53;
posx_h=0.648;
posz_h=1.3;
%y_afil=0.55;
posx_gmp=0.567;
posz_gmp=0;
airfoil_number=4;

%Monoplano V1
% Tr1_w=0.8137;
% Tr2_w=0.7469;
% b_h=1.55;
% b_w=2.3;
% cr_h=0.344;
% cr_w=0.51;
% posx_h=0.52;
% posz_h=0.32;
% y_afil=0.56;
% posx_gmp=0.51;
% posz_gmp=0.225;
% airfoil_number=1;

%Wing(s)
iswing=true;

geo.nwing=1;                    %number of wings (scalar)
geo.nelem=[2];				%number of partitions on each wing (1d array)
geo.ref_point=[0 0 0];      %System reference point
%geo.fnx(s,t)=[~,~];
geo.c=[cr_w];					%Root chord (2d array) [wing1 wing2]
geo.symetric=[1];			    %Wing symmetry boolean bit (2d array) [wing1 wing2]
geo.startx=[geo.ref_point(1)-geo.c/4];		    	%Partition starting coordinate (2d array) [wing1 wing2]
geo.starty=[0];		    	% ---"----
geo.startz=[0];		    	% ---"----

% geo.c=[1 0.5];					%Root chord (2d array) [wing1 wing2]
% geo.symetric=[1 1];			    %Wing symmetry boolean bit (2d array) [wing1 wing2]
% geo.startx=[0 2];		    	%Partition starting coordinate (2d array) [wing1 wing2]
% geo.starty=[0 0];		    	% ---"----
% geo.startz=[0 1];		    	% ---"----

% Wing 1
geo.T(1,:)=[Tr1_w Tr1_w*Tr2_w];					%Taper ratio (2d array)
geo.b(1,:)=[0.55*b_w/2 (1-0.55)*b_w/2];					%span(distance root->tip chord) (2d array)
geo.SW(1,:)=[-atan(geo.c*(1-geo.T(1))/(4*geo.b(1))) -atan(geo.c*(geo.T(1)-geo.T(2))/(4*geo.b(2)))];					%Sweep (2d array)
geo.dihed(1,:)=[0 0];				%Dihedral (2d array)	
geo.fnx(1,:)=[0 0];                 %number of panels on flap chords (2d array)
geo.ny(1,:)=[5 5];					%number of panels in span (2d array)
geo.nx(1,:)=[5 5];					%number of panels on chord (2d array)
geo.fsym(1,:)=[0 0];				%flap deflection symmetry boolean bit  (2d array)
geo.fc(1,:)=[0 0];					%flap chord in percent of wingchord (2d array)
geo.flapped(1,:)=[0 0];			    %flapped partition(wing part) boolean bit (2d array)
geo.TW(1,1:2,1)=0;                  %partition twist (3d array)<1 inboard, 2 outboard>
geo.TW(1,1:2,2)=0;

Perfil_af1=char(airfoil_lib(airfoil_number));
Perfil_af2=char(airfoil_lib(airfoil_number));

geo.foil(1,1:2,1)={strcat(Perfil_af1,'.dat')};		    %Partition airfoils (3d array)	
geo.foil(1,1,2)={strcat(Perfil_af1,'.dat')};              %1 inboard, 2 outboard
geo.foil(1,2,2)={strcat(Perfil_af2,'.dat')};    %1 inboard, 2 outboard

% Wing 2
% geo.T(2,:)=[1 0];					%Taper ratio (2d array)
% geo.b(2,:)=[2 0];					%span(distance root->tip chord) (2d array)
% geo.SW(2,:)=[0 0];					%Sweep (2d array)
% geo.dihed(2,:)=[0 0];				%Dihedral (2d array)	
% geo.fnx(2,:)=[0 0];
% geo.ny(2,:)=[5 0];					%number of panels in span (2d array)
% geo.nx(2,:)=[5 0];					%number of panels on chord (2d array)
% geo.fsym(2,:)=[0 0];				    %flap deflection symmetry boolean bit  (2d array)
% geo.fc(2,:)=[0 0];					%flap chord in percent of wingchord (2d array)
% geo.flapped(2,:)=[0 0];			    %flapped partition(wing part) boolean bit (2d array)
% geo.TW(2,:,1)=0;	    	%partition twist (3d array)<1 inboard, 2 outboard>
% geo.TW(2,:,2)=0;
% geo.foil(2,:,1)={'2412'};		    %Partition airfoils (3d array)	
% geo.foil(2,:,2)={'2412'};		    %1 inboard, 2 outboard

geo.flap_vector=zeros(size(geo.flapped));          %Flap deflection vector 
geo.CG=[0 0 0];             %System center of gravity (around which all rotations occur)
%geo.CG=[0.009 0 0.1];             %System center of gravity (around which all rotations occur)

state.AS=15;					%airspeed
state.alpha=0;				%angle of attack
state.betha=0;				%angle of sideslip
state.P=0;					%roll angluar rate	
state.Q=0;					%pitch angular rate
state.R=0;					%Yaw angular rate
state.alphadot=0;           %Angle of attack time derivative
state.bethadot=0;           %Angle of sidesliptime derivative
state.ALT=0;                %Altitude, meters.
state.rho=1.1;                %Air density, kg/m^3.
state.pgcorr=0;             %Prandtl-Glauert compressibillity correction.

lattice.XYZ=0;				%panel corner matrix (2d array)
lattice.COLLOC=0;           %collocation point matrix
lattice.VORTEX=0;           %Vortex sling cornerpoint position matrix
lattice.N=0;                %Airfoil collocation point normal direction matrix

ref.S_ref=0;                %reference area;
ref.C_mac=0;                %mean aerodynamic choord
ref.mac_pos=0;              %start position of mac
ref.C_mgc=0;                %mean geometric chord
ref.b_ref=0;                %reference span

results.dwcond=0;           %computation result memory structure.


Polar_aux1=load(strcat('D:\AeroDesign\ALPHA 2020\MDO 2020\Validação\Polares\Polar_',Perfil_af1,'.txt'));
Polar_aux2=load(strcat('D:\AeroDesign\ALPHA 2020\MDO 2020\Validação\Polares\Polar_',Perfil_af2,'.txt'));

if max(Polar_aux1(:,2))>max(Polar_aux2(:,2))
    Polar2D=strcat('Polar_',Perfil_af1);
else
    Polar2D=strcat('Polar_',Perfil_af2);
end
    
[wing.polar,ref]=Aero3D(geo,state,lattice,ref,results,alpha_inicial,alpha_incremento,alpha_final,Polar2D,iswing);

wing.geo=geo;
wing.ref=ref;
wing.state=state;

% fPath = 'D:\Faculdade\T135_export\Results\';
% fName = strcat('wing.mat');
% save(strcat(fPath, fName))

%--------------------Horizontal Stabilizer-----------------------------
clear geo
clear ref
clear lattice

%modeFrontier Stabilizer
geo.b(1,:)=[b_h/2];					%span(distance root->tip chord) (2d array)
geo.c=[cr_h];					%Root chord (2d array) [wing1 wing2]
geo.ref_point=[posx_h 0 posz_h]; %System reference point

geo.startx=[geo.ref_point(1)-geo.c/4];		    	%Partition starting coordinate (2d array) [wing1 wing2]
geo.startz=[geo.ref_point(3)];		    	% ---"----
geo.starty=[0];		    	% ---"----



alpha_inicial= min(wing.polar{1,1}(:,1))*180/pi;
alpha_incremento=1;
alpha_final=max(wing.polar{1,1}(:,1))*180/pi;

iswing=false;

geo.nwing=1;                    %number of wings (scalar)
geo.nelem=[1];				%number of partitions on each wing (1d array)
geo.symetric=[1];			    %Wing symmetry boolean bit (2d array) [wing1 wing2]



geo.T(1,:)=[1];					%Taper ratio (2d array)
geo.SW(1,:)=[0];					%Sweep (2d array)
geo.dihed(1,:)=[0];				%Dihedral (2d array)	
geo.fnx(1,:)=[0];                 %number of panels on flap chords (2d array)
geo.ny(1,:)=[5];					%number of panels in span (2d array)
geo.nx(1,:)=[5];					%number of panels on chord (2d array)
geo.fsym(1,:)=[0];				%flap deflection symmetry boolean bit  (2d array)
geo.fc(1,:)=[0];					%flap chord in percent of wingchord (2d array)
geo.flapped(1,:)=[0];			    %flapped partition(wing part) boolean bit (2d array)
geo.TW(1,:,1)=0;                  %partition twist (3d array)<1 inboard, 2 outboard>
geo.TW(1,:,2)=0;
geo.flap_vector=zeros(size(geo.flapped));          %Flap deflection vector 



geo.CG=wing.geo.CG;             %System center of gravity (around which all rotations occur)

lattice.XYZ=0;				%panel corner matrix (2d array)
lattice.COLLOC=0;           %collocation point matrix
lattice.VORTEX=0;           %Vortex sling cornerpoint position matrix
lattice.N=0;                %Airfoil collocation point normal direction matrix

ref.S_ref=0;                %reference area;
ref.C_mac=0;                %mean aerodynamic choord
ref.mac_pos=0;              %start position of mac
ref.C_mgc=0;                %mean geometric chord
ref.b_ref=0;                %reference span

results.dwcond=0;           %computation result memory structure.



foil_stabilizer='1598_'; %nome do perfil
j=1;
delta_deflex=5;

for i=-30:delta_deflex:30
   geo.foil(1,:,1)={strcat(foil_stabilizer,num2str(i),'.dat')};		    %Partition airfoils (3d array)	
   geo.foil(1,:,2)={strcat(foil_stabilizer,num2str(i),'.dat')};         %1 inboard, 2 outboard
   Polar2D=strcat('Polar_1598_',num2str(i));
   [stabilizer.polar(j),ref,k]=Aero3D(geo,state,lattice,ref,results,alpha_inicial,alpha_incremento,alpha_final,Polar2D,iswing);
   j=j+1;
end   
stabilizer.n_deflex=j-1; %nº de deflexões 
stabilizer.geo=geo;
stabilizer.ref=ref;

% fPath = 'D:\Faculdade\T135_export\Results\';
% fName = strcat('stabilizer.mat');
% save(strcat(fPath, fName))

% -------------------- GMP ----------------------------------------------

GMP.T_stat=43.377;  % Tração estática [N]
GMP.d_t1=-0.6843; %Derivada da tração de 1º grau [-]
%GMP.d_t2=-0.0242; %Derivada da tração de 2º grau [-]
GMP.jp=0.19481;       %Inércia da hélice
GMP.diam_p=0.35; %diametro da hélice [m]
GMP.n_hel=1;    %nº de pás da hélice [-]
GMP.n_motor=1;    %nº de motores [-]
GMP.posx=posx_gmp;
GMP.ip=0;       %Incidencia do motor [º]
GMP.posz=posz_gmp;
GMP.rot=13200;   %RPM
%[Cm_aircraft,Cd_aircraft,Cl_aircraft,epslon,Max_stab,ME_min,Cl_trim,Cd_trim] = Stability(wing,stabilizer,GMP);
%[aircraft,downwash,Max_stab,ME_min,Cl_trim,Cd_trim,alpha_limite] = Stability( wing,stabilizer,GMP );
[ Matrix,results ] =  Ideal_CG( wing,stabilizer,GMP,xCG,zCG,constrains);

MTOW=results.MTOW;
AR=(wing.ref.b_ref^2)/wing.ref.S_ref;
S_w=wing.ref.S_ref;
S_h=stabilizer.ref.S_ref;
CG_x=results.CG(1);
CG_z=results.CG(3);