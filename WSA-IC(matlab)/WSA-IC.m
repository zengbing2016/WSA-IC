%{
    WSA-IC for multimodal function optimization
    Liang Gao(email: gaoliang@mail.hust.edu.cn)
    Bing Zeng(email: zengbing2016@126.com)
%}

clear;clc;

n=5;iter=300000;
Ts=100*n;Tf=1.0*10^-8;
size=60;
vLBd=-100;vUBd=100;
fgbest=inf;
GloOpt=[];GloFit=[];

for i=1:1:size
    E.value(i,:)=rand(1,n)*(vUBd-vLBd)+vLBd;
    E.c(i)=0;
    E.fitness(i)=Func(E.value(i,:));
    if (E.fitness(i)<fgbest)
        v=i;
        fgbest=E.fitness(i);
    end
end

m=0;
while (m<iter)
    m=m+1
    for i=1:1:size
        v=FindBetterNearest(i, E(:,:),size);
        if (v~=-1)
            X=E.value(i,:)+rand(1,n)*2.*(E.value(v,:)-E.value(i,:));
            X((X<vLBd))=vLBd;
            X((X>vUBd))=vUBd;
            T=Func(X);
            if (E.fitness(i)>T)
                E.value(i,:)=X(:);
                E.fitness(i)=T;
                E.c(i)=0;
            else
                if (E.c(i)~=Ts)
                    E.c(i)=E.c(i)+1;
                else
                    [GloOpt,GloFit,fgbest] = JudWthGloOpt(GloOpt,GloFit,fgbest,Tf,E,i);
                    E.value(i,:)=rand(1,n)*(vUBd-vLBd)+vLBd;
                    E.c(i)=0;
                    E.fitness(i)=Func(E.value(i,:));
                end
            end
        else
            if (E.c(i)~=Ts)
                E.c(i)=E.c(i)+1;
            else
                [GloOpt,GloFit,fgbest] = JudWthGloOpt(GloOpt,GloFit,fgbest,Tf,E,i);
                E.value(i,:)=rand(1,n)*(vUBd-vLBd)+vLBd;
                E.c(i)=0;
                E.fitness(i)=Func(E.value(i,:));
            end
        end
    end
end

for i=1:1:size
    [GloOpt,GloFit,fgbest] = JudWthGloOpt(GloOpt,GloFit,fgbest,Tf,E,i);
end

Td=28;
GloOptFinal=[];GloFitFinal=[];
for i=1:1:length(GloFit)
    if abs(GloFit(i)-fgbest)<Tf
        Judge=0;
        for j=1:1:length(GloFitFinal)
            if norm((GloOptFinal(j,:)-GloOpt(i,:)),2)<Td
                Judge=1;
                break;
            end
        end
        if (Judge==0)
            row=length(GloFitFinal)+1;
            GloOptFinal(row,:)=GloOpt(i,:);
            GloFitFinal(row)=GloFit(i);
        end
    end
end

name='global_optima';
save(name,'GloOptFinal','GloFitFinal');