function output=FindBetterNearest(n,E,size)
Min=inf;output=-1;
for i=1:1:size
    if (E.fitness(i)<E.fitness(n))
        M=norm((E.value(i,:)-E.value(n,:)),2);
        if (M<Min)
            Min=M;
            output=i;
        end
    end
end
end