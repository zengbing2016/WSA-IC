function ObjVal=Func(X)
%Five-Uneven-Peak Trap
ObjVal=0;
for i=1:1:length(X)
    if (X(i)<0)
        ObjVal=ObjVal-200.0+X(i)^2;
    else if (X(i)<2.5)
            ObjVal=ObjVal-80.0*(2.5-X(i));
        else if (X(i)<5.0)
                ObjVal=ObjVal-64.0*(X(i)-2.5);
            else if (X(i)<7.5)
                    ObjVal=ObjVal-160.0+X(i)^2;
                else if (X(i)<12.5)
                        ObjVal=ObjVal-28.0*(X(i)-7.5);
                    else if (X(i)<17.5)
                            ObjVal=ObjVal-28.0*(17.5-X(i));
                        else if (X(i)<22.5)
                                ObjVal=ObjVal-32.0*(X(i)-17.5);
                            else if (X(i)<27.5)
                                    ObjVal=ObjVal-32.0*(27.5-X(i));
                                else if (X(i) <= 30.0)
                                        ObjVal=ObjVal-80.0*(X(i)-27.5);
                                    else
                                        ObjVal=ObjVal-200.0+(X(i)-30.0)^2;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

ObjVal=ObjVal+200.0*length(X);

%rosenbrock
% dim=length(X);
% mat0=X(:,1:dim-1);
% mat1=X(:,2:dim);
% if dim == 2
% 	ObjVal=100*(mat1-mat0.^2).^2+(1-mat0).^2;
% else
% 	ObjVal=sum(100*(mat1-mat0.^2).^2+(1-mat0).^2,2);
% end

return