function [ polar_corrigida,downwash,alpha_0w ] = Polar_correction( wing,stabilizer )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

z_h=stabilizer.geo.ref_point(3);
z_w=wing.geo.ref_point(3);
x_h=stabilizer.geo.ref_point(1);
x_w=wing.geo.ref_point(1);
b_w=wing.ref.b_ref;
Cr_w=wing.geo.c;
Ct_w=Cr_w*wing.geo.T(1,2);
S_w=wing.ref.S_ref;
Sweep_w=wing.geo.SW(1,1);
alpha_w=wing.polar{1,1}(:,1);
AR_w=(b_w^2)/S_w;

%% Alpha para zero lift
in=find(~isnan(wing.polar{1,1}(:,2)), 1);
out=find(~isnan(wing.polar{1,1}(:,2)), 1, 'last');
alpha_0w=interp1(wing.polar{1,1}(in:out,2),wing.polar{1,1}(in:out,1),0,'pchip');

%% Downwash %top
Kh=(1-((z_h-z_w)/b_w))/(2*((x_h-x_w)/b_w))^(1/3);
Kdelta=(10-3*Ct_w/Cr_w)/7;
Ka=(1/AR_w)-1/(1+AR_w^1.7);

deda=4.44*(Ka*Kdelta*Kh*sqrt(cos(deg2rad(Sweep_w))))^1.19;
downwash=deda*(alpha_w-alpha_0w);

alpha_ef=alpha_w-downwash;

for j=1:stabilizer.n_deflex
    [row,~]=find(~isnan(stabilizer.polar{1,j}(:,2)));
    pos_min=min(row);
    pos_max=max(row);
    stabilizer.polar{1,j}(pos_min:pos_max,2)=interp1(stabilizer.polar{1,j}(pos_min:pos_max,1),stabilizer.polar{1,j}(pos_min:pos_max,2),alpha_ef(pos_min:pos_max),'pchip'); %CL
    stabilizer.polar{1,j}(pos_min:pos_max,3)=interp1(stabilizer.polar{1,j}(pos_min:pos_max,1),stabilizer.polar{1,j}(pos_min:pos_max,3),alpha_ef(pos_min:pos_max),'pchip'); %CDi
    stabilizer.polar{1,j}(pos_min:pos_max,4)=interp1(stabilizer.polar{1,j}(pos_min:pos_max,1),stabilizer.polar{1,j}(pos_min:pos_max,4),alpha_ef(pos_min:pos_max),'pchip'); %CDp
    stabilizer.polar{1,j}(pos_min:pos_max,5)=interp1(stabilizer.polar{1,j}(pos_min:pos_max,1),stabilizer.polar{1,j}(pos_min:pos_max,5),alpha_ef(pos_min:pos_max),'pchip'); %CM 
%    stabilizer.polar{1,j}(:,6)=alpha_ef(pos_min:pos_max)*180/pi;
end
polar_corrigida=stabilizer.polar;
end

