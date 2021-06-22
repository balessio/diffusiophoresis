% solves for y-averaged colloid concentration using 1D approximate equation
global time_now
global B c0 l_h Dp_Ds 
B = 0.1; c0 = 1.0; l_h = 10; Dp_Ds = 1e-2;

Peclet = 0.334;
PsiW = -4; PsiD = -3;
D_pos_D_Na = 1.0; D_neg_D_Na = 1.526;
global Gw_Ds Gp_Ds
Gw_Ds = Peclet*(((D_pos_D_Na-D_neg_D_Na)/(D_pos_D_Na+D_neg_D_Na))*PsiW+PsiW^2/8);
Gp_Ds = Peclet*(((D_pos_D_Na-D_neg_D_Na)/(D_pos_D_Na+D_neg_D_Na))*PsiD+PsiD^2/8);

start_time = 0.001;
global erf_denom
erf_denom = sqrt(4*start_time);
if start_time==0.0
   erf_denom = 0.01; 
end

%x = [linspace(0,0.09,100) linspace(0.1,1,10)];
x = linspace(0,1,50);
t = linspace(start_time,1.05,25);

[cSolute, dcSolutedX, d2cSolutedX2] = solute_withDerivatives(x,time_now);

VDP = Gp_Ds.*dcSolutedX./cSolute; Vs = Gw_Ds.*dcSolutedX./cSolute;
dVDPdX = Gp_Ds*((cSolute.*d2cSolutedX2 - dcSolutedX.^2)./cSolute.^2);

calDp = 1 + (2/105)*l_h^-2*Dp_Ds^-2.*Vs.*Vs;



figure
set(gcf,'Position',[800 800 800 800])
suptitle(['time = ' num2str(time_now)])
subplot(4,2,1)
plot(x,cSolute)
ylabel('c')
subplot(4,2,2)
plot(x,dcSolutedX)
ylabel('dc/dx')
subplot(4,2,3)
plot(x,d2cSolutedX2)
ylabel('d2c/dx2')
subplot(4,2,4)
plot(x,VDP)
ylabel('V_{DP}')
subplot(4,2,5)
plot(x,dVDPdX)
ylabel('dV_{DP}/dx')
subplot(4,2,6)
plot(x,Vs)
ylabel('V_s')
subplot(4,2,7)
plot(x,calDp)
ylabel('calD_p')


function [C,dCdx,d2Cdx2] = solute_withDerivatives(xx,tt)
global B

res=1.0;
C=zeros(size(xx));
dCdx=zeros(size(xx));
d2Cdx2=zeros(size(xx));

i=0;
while(res>1.0e-4)
lambda=(2*i+1)*pi/2;
val=C;
C = C + (1-B)*2*(1-cos(lambda))/lambda*sin(lambda*xx)*exp(-lambda^2*tt);
res=sqrt(sum(sum((C-val).^2))./sum(sum(C.^2)));
dCdx = dCdx + (1-B)*2*(1-cos(lambda))*cos(lambda*xx)*exp(-lambda^2*tt);
d2Cdx2 = d2Cdx2 -(1-B)*2*(1-cos(lambda))*lambda*sin(lambda*xx)*exp(-lambda^2*tt);
i=i+1;
end
C=C+B;
end