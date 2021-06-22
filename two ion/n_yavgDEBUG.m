% solves for y-averaged colloid concentration using 1D approximate equation

global B c0 l_h Dp_Ds 
B = 0.6; c0 = 1.0; l_h = 10; %Dp_Ds = 1e-2;

start_time = 0.0;

Peclet = 0.334;
PsiW = -4; PsiD = -3;
%global PsiW
%PsiD = -3;
D_pos_D_Na = 1.0; D_neg_D_Na = 1.526;
global Gw_Ds Gp_Ds
Gw_Ds = Peclet*(((D_pos_D_Na-D_neg_D_Na)/(D_pos_D_Na+D_neg_D_Na))*PsiW+PsiW^2/8);
Gp_Ds = Peclet*(((D_pos_D_Na-D_neg_D_Na)/(D_pos_D_Na+D_neg_D_Na))*PsiD+PsiD^2/8);


global erf_denom
erf_denom = sqrt(4*start_time);
if start_time==0.0
   erf_denom = 0.01; 
end

%x = [linspace(0,0.09,100) linspace(0.1,1,10)];
x = linspace(0,1,25);
t = linspace(start_time,1,25);

A1 = (1/105)*(1/(l_h*Dp_Ds)^2)*(7*Gp_Ds*Gw_Ds-4*Gw_Ds^2);
A2 = (1/105)*(1/(l_h*Dp_Ds)^2)*(7*Gp_Ds*Gw_Ds-2*Gw_Ds^2);

global time_now

[cSolute, dcSolutedX, d2cSolutedX2] = solute_withDerivatives(x,time_now);
Up_Ds = Gp_Ds.*dcSolutedX./cSolute + A1.*(dcSolutedX./cSolute.^3).*(cSolute.*d2cSolutedX2-dcSolutedX.^2);
calDp_Ds = Dp_Ds - A2.*(dcSolutedX./cSolute).^2;

c = calDp_Ds.^-1;
s = -Up_Ds./calDp_Ds;


%plot(x,cSolute,'linewidth',3)
%plot(x,dcSolutedX,'linewidth',3)
%plot(x,d2cSolutedX2,'linewidth',3)
%plot(x,c,'linewidth',3)
%plot(x,s,'linewidth',3)
plot(x,Up_Ds,'linewidth',3)
%plot(x,calDp_Ds,'linewidth',3)


%----------------------------------------------

function [C,dCdx,d2Cdx2] = solute_withDerivatives(xx,tt)
global B

res=1.0;
C=zeros(size(xx));
dCdx=zeros(size(xx));
d2Cdx2=zeros(size(xx));

i=0;
while(res>1.0e-2)
lambda=(2*i+1)*pi/2;
blambda = 2*(1-cos(lambda));
val=C;
C = C + (1-B)*blambda/lambda*sin(lambda*xx)*exp(-lambda^2*tt);
res=sqrt(sum(sum((C-val).^2))./sum(sum(C.^2)));
dCdx = dCdx + (1-B)*blambda*cos(lambda*xx)*exp(-lambda^2*tt);
d2Cdx2 = d2Cdx2 -(1-B)*blambda*lambda*sin(lambda*xx)*exp(-lambda^2*tt);
i=i+1;
end
C=C+B;
end