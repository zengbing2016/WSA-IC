function [GloOpt,GloFit,fgbest] = JudWthGloOpt(GloOpt,GloFit,fgbest,Tf,E,i)
if (E.fitness(i)<fgbest)
    if (abs(fgbest-E.fitness(i))>Tf)
        GloOpt=[];
        GloFit=[];
    end
    fgbest=E.fitness(i);
    row=size(GloOpt,1)+1;
    GloOpt(row,:)=E.value(i,:);
    GloFit(row)=E.fitness(i);
else
    if(abs(E.fitness(i)-fgbest)<=Tf)
        row=size(GloOpt,1)+1;
        GloOpt(row,:)=E.value(i,:);
        GloFit(row)=E.fitness(i);
    end
end
end