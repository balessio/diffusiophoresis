function [C,dCdx] = solute_withDerivative_copy(xx,tt)
global B

res=1.0;
C=zeros(size(xx));
dCdx=zeros(size(xx));

i=0;
while(res>1.0e-4)
lambda=(2*i+1)*pi/2;
val=C;
C = C + (1-B)*2*(1-cos(lambda))/lambda*sin(lambda*xx)*exp(-lambda^2*tt);
res=sqrt(sum(sum((C-val).^2))./sum(sum(C.^2)));
dCdx = dCdx + (1-B)*2*(1-cos(lambda))*cos(lambda*xx)*exp(-lambda^2*tt);
i=i+1;
end
C=C+B;
end