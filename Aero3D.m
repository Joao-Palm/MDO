function [Polar,ref,j] = Aero3D(geo,state,lattice,ref,results,alpha_inicial,alpha_incremento,alpha_final,Polar2D,iswing)


Polar2D=load(strcat('D:\AeroDesign\ALPHA 2020\MDO 2020\Validação\Polares\',Polar2D,'.txt'));
results=[];
[CLmax,idmax] = max(Polar2D(:,2));
[CLmin, idmin] = min(Polar2D(:,2));

a1=alpha_inicial*pi/180;
b1=alpha_incremento*pi/180;
c1=alpha_final*pi/180;
if ~iswing
    temp=a1;
    a1=c1;
    c1=temp;
    b1=-b1;
end
j=0;
[lattice,ref]=fLattice_setup2(geo,state,1);


for alpha=a1:b1:c1
    
    state.alpha=alpha;
    j=j+1;
    
    
    [results]=solver9(results,state,geo,lattice,ref);
    [results]=coeff_create3(results,lattice,state,ref,geo);
    
    
    %results.alpha_sweep(j)=state.alpha;
%     if alpha<0
%         CDptemp = interp1(Polar2D(idmin:idmax, 2), Polar2D(idmin:idmax, 4), results.CLwing(1),'linear','extrap');%+interp1(Polar2D(idmin:idmax, 2), Polar2D(idmin:idmax, 4), results.CLwing(2));    
%     else
%         CDptemp = interp1(Polar2D(idmin:idmax, 2), Polar2D(idmin:idmax, 4), results.CLwing(1));    
%     end
    CDptemp = interp1(Polar2D(idmin:idmax, 2), Polar2D(idmin:idmax, 4), results.CLwing(1));    

%    CDptemp = interp1(Polar2D(idmin:idmax, 2), Polar2D(idmin:idmax, 4), results.CLwing(1));
    results.Polar(j,:)=[state.alpha results.CL results.CD CDptemp results.Cm results.CL_a];
    %results.Polar(j,:)=[state.alpha results.CL results.CD CDptemp results.Cm results.CL_a];
    
    if iswing
        if max(max(results.CL_local)) > CLmax
            stall = true;
            break
        end
        
        if max(results.CLwing) > CLmax
            stall = true;
            break;
        end
        % Se o CL é menor do que pro angulo anterior
        if j > 1 && results.CL <  results.Polar(j-1,2)
            stall = true;
            break;
        end
    else
%                      if min(min(results.CL_local)) < CLmin
%                          stall = true;
%                          break
%                      end
        
%         if (results.CLwing) < CLmin
%             stall = true;
%             break;
%         end
    end
    
    
    if ~isnan(CDptemp)
        results.Polar(j,4) = CDptemp;
    else
        results.Polar(j,:) = NaN;
        %break
    end%CDp
    
end

if ~iswing
    matriz_aux=results.Polar;
    for i=1:length(results.Polar(:,1))
        results.Polar(i,:)=matriz_aux(length(results.Polar(:,1))+1-i,:);
    end
end

%postproc(lattice,geo,ref)
Polar={results.Polar};


end

