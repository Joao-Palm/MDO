function [ MTOW,VDec,CLdec] = Performance2( wing,stabilizer,GMP,Cl_trim,downwash,alpha_limite)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

alpha=wing.polar{1,1}(1:alpha_limite,1);

m = 8:0.2:50;
ax_min = 0;
decLim = 50;
g=9.81;
rho=1.2;
mi_rol=0.065;

[CLdec, indexDec] = max(Cl_trim(:,1));
in=find(~isnan(wing.polar{1,1}(:,2)), 1);
out=find(~isnan(wing.polar{1,1}(:,2)), 1, 'last');

AlphaDec = alpha(indexDec,1);
AlphaDec_w = AlphaDec + wing.geo.TW(1,1,1);
AlphaDec_h = AlphaDec-interp1(alpha(in:alpha_limite,1),downwash(in:alpha_limite,1),AlphaDec_w,'pchip');% + AlphaDec_w;

AlphaCor_w = wing.geo.TW(1,1,1);
AlphaCor_h = -interp1(alpha(in:alpha_limite,1),downwash(in:alpha_limite,1),AlphaCor_w,'pchip');

take_off = true;
if CLdec < 0
    take_off = false;
end

%Wing Polar
CLdec_w=interp1(wing.polar{1,1}(in:out, 1),wing.polar{1,1}(in:out, 2),AlphaDec_w, 'pchip');
CDdec_w=interp1(wing.polar{1,1}(in:out, 1),wing.polar{1,1}(in:out, 3),AlphaDec_w, 'pchip')+interp1(wing.polar{1,1}(in:out, 1),wing.polar{1,1}(in:out, 4),AlphaDec_w, 'pchip');
CMdec_w=interp1(wing.polar{1,1}(in:out, 1),wing.polar{1,1}(in:out, 4),AlphaDec_w, 'pchip');

CLcor_w=interp1(wing.polar{1,1}(in:out, 1),wing.polar{1,1}(in:out, 2),AlphaCor_w, 'pchip');
CDcor_w=interp1(wing.polar{1,1}(in:out, 1),wing.polar{1,1}(in:out, 3),AlphaCor_w, 'pchip');
CMcor_w=interp1(wing.polar{1,1}(in:out, 1),wing.polar{1,1}(in:out, 4),AlphaCor_w, 'pchip');

%Stabilizer Polar
[row,~]=find(~isnan(stabilizer.polar{1,7}(:,2)));
pos_min=min(row);
pos_max=max(row);
CLdec_h=interp1(stabilizer.polar{1,7}(pos_min:pos_max,1),stabilizer.polar{1,7}(pos_min:pos_max,2),AlphaDec_h, 'pchip');
CDdec_h=interp1(stabilizer.polar{1,7}(pos_min:pos_max,1),stabilizer.polar{1,7}(pos_min:pos_max,3),AlphaDec_h, 'pchip')+interp1(stabilizer.polar{1,7}(pos_min:pos_max,1),stabilizer.polar{1,7}(pos_min:pos_max,4),AlphaDec_h, 'pchip');
CMdec_h=interp1(stabilizer.polar{1,7}(pos_min:pos_max,1),stabilizer.polar{1,7}(pos_min:pos_max,5),AlphaDec_h, 'pchip');

CLcor_h=interp1(stabilizer.polar{1,7}(pos_min:pos_max,1),stabilizer.polar{1,7}(pos_min:pos_max,2),AlphaCor_h, 'pchip');
CDcor_h=interp1(stabilizer.polar{1,7}(pos_min:pos_max,1),stabilizer.polar{1,7}(pos_min:pos_max,3),AlphaCor_h, 'pchip')+interp1(stabilizer.polar{1,7}(pos_min:pos_max,1),stabilizer.polar{1,7}(pos_min:pos_max,4),AlphaCor_h, 'pchip');
CMcor_h=interp1(stabilizer.polar{1,7}(pos_min:pos_max,1),stabilizer.polar{1,7}(pos_min:pos_max,5),AlphaCor_h, 'pchip');

if take_off == true
 %           V_dec = sqrt((2*m*g)/(rho * (wing.ref.S_ref *CLdec_w + stabilizer.ref.S_ref *CLdec_h)));

            V_dec = sqrt((2*m*g)./(rho * (wing.ref.S_ref *CLdec)));
            dec.a = -0.5*rho./m * (wing.ref.S_ref*(CDcor_w - mi_rol * CLcor_w) + ...
                stabilizer.ref.S_ref*(CDcor_h - mi_rol * CLcor_h));

            dec.b = GMP.d_t1 * (cos(GMP.ip*pi/180) + mi_rol *sin(GMP.ip*pi/180)) ./m;

            dec.c = -1./m .* (mi_rol*g .* m - GMP.T_stat*(cos(GMP.ip*pi/180) + mi_rol*sin(GMP.ip*pi/180)));
            
            dec.Delta = dec.b.^2 - 4*dec.a .* dec.c;
            deltaNeg = find(dec.Delta < 0);
            dec.Delta(deltaNeg) = [];
            dec.a(deltaNeg) = [];
            dec.b(deltaNeg) = [];
            dec.c(deltaNeg) = [];
            V_dec(deltaNeg) = [];
            m(deltaNeg) = [];
            
            dec.r1 = (-dec.b + sqrt(dec.Delta))./(2*dec.a);
            dec.r2 = (-dec.b - sqrt(dec.Delta))./(2*dec.a);
            dec.U1 = dec.r1./sqrt(dec.Delta);
            dec.U2 = -dec.r2./sqrt(dec.Delta);
            
            if isreal(dec.r1)
                dec.S  = dec.U2 .* log(abs(1 - V_dec./dec.r2)) + dec.U1 .* log(abs(1 - V_dec./dec.r1));
            else
                dec.S = -1;
            end
           
            decCorte = find(dec.S > decLim, 1, 'first');
            dec.S(decCorte+1:end) = [];
            m(decCorte+1:end) = [];
            V_dec(decCorte+1:end) = [];
            q_dec = 0.5 * rho * V_dec.^2;

            if length(dec.S) < 2
                MTOWcl = 0;
                MTOWax = 0;
            else
                % Se nunca decola por CL
                if max(dec.S) < 0 || dec.S(1)-dec.S(2) > 0
                    MTOWcl = 0;
                    MTOWax = 0;
                    % Se decola por CL entre 0 e 50m
                else
                    MTOWcl = interp1(dec.S, m, decLim);
                    dec.Fx_w   = -q_dec * wing.ref.S_ref * CDdec_w;
                    dec.Fx_h   = -q_dec * stabilizer.ref.S_ref * CDdec_h;
                    dec.Fx_gmp = (GMP.T_stat + GMP.d_t1 * V_dec) * cos(GMP.ip*pi/180);
                    dec.Fx = dec.Fx_w + dec.Fx_h + dec.Fx_gmp;
                    dec.ax = dec.Fx./m;
                        % Se aceleracao máxima é menor do que o minimo imposto
                    if max(dec.ax) < ax_min
                        MTOWax  = 0;
                        % Se aceleracao minima é maior do que o minimo imposto
                    elseif min(dec.ax) > ax_min
                        MTOWax = max(m);
                        % C.C.
                    else
                        MTOWax = interp1(dec.ax, m, ax_min);
                    end
                end
            end

            MTOW = min(MTOWcl, MTOWax);
end
    
VDec=sqrt((2*MTOW*g)/(rho * (wing.ref.S_ref *CLdec_w + stabilizer.ref.S_ref *CLdec_h)));
% a=wing.ref.S_ref *CLdec_w + stabilizer.ref.S_ref *CLdec_h
% b=wing.ref.S_ref *CLdec
end

