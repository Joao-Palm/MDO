function [ Matrix,results ] =  Ideal_CG( wing,stabilizer,GMP,xCG,zCG,constrains)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[ stabilizer.polar,downwash,alpha_0w ] = Polar_correction( wing,stabilizer );
for i=1:length(xCG)
    for j=1:length(zCG)
        wing.geo.CG=[xCG(i) 0 zCG(j)];
        [Max_stab,ME_min,Cl_trim,alpha_limite]=Stability(wing,stabilizer,GMP,downwash,alpha_0w);
        if ~isnan(max(Cl_trim(:,1)))
            [ MTOW,VDec,CLdec] = Performance2( wing,stabilizer,GMP,Cl_trim,downwash,alpha_limite);
        else
            MTOW=NaN;
            VDec=NaN;
            CLdec=NaN;
        end
        Matrix.MTOW(i,j)=MTOW;
        Matrix.ME_min(i,j)=ME_min;
        Matrix.Max_stab(i,j)=Max_stab;
        Matrix.VDec(i,j)=VDec;
        Matrix.Cltrim(i,j)=CLdec;
    end
end


ME_filtrada=Matrix.ME_min;
ME_filtrada(ME_filtrada<constrains.ME)=NaN;

Max_stab_filtrada=Matrix.Max_stab;
Max_stab_filtrada(Max_stab_filtrada>constrains.Stab)=NaN;

for i=1:length(xCG)
    for j=1:length(zCG)
        if (Matrix.ME_min(i,j)<constrains.ME) || (Matrix.Max_stab(i,j)>constrains.Stab)
            MTOW_filtrada(i,j)=NaN;
        else
            MTOW_filtrada(i,j)=Matrix.MTOW(i,j);
        end
    end
end

[results.MTOW,pos]=max(MTOW_filtrada(:));
[results.pos(1), results.pos(2)] = ind2sub(size(MTOW_filtrada),pos);
results.CG=[xCG(results.pos(1)) 0 zCG(results.pos(2))];
results.ME=ME_filtrada(results.pos(1),results.pos(2));
results.Max_stab=Max_stab_filtrada(results.pos(1),results.pos(2));
results.VDec=Matrix.VDec(results.pos(1),results.pos(2));
results.Cl_trim=Matrix.Cltrim(results.pos(1),results.pos(2));
end

