% solves for y-averaged colloid concentration using 1D approximate equation

global B c0 l_h Dp_Ds 
B = 0.1; c0 = 1.0; l_h = 10; Dp_Ds = 1e-2;

Peclet = 0.334;
%global PsiW
PsiW = -4;
PsiD = -3;
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
x = linspace(0,1,51);
t = linspace(start_time,2.05,205);

m = 0;
N_yavg_solved = pdepe(m,@N_yavg,@initial_conditions,@boundary_conditions,x,t);

surf(x,t,N_yavg_solved(:,:,1))
xlabel('x')
ylabel('t')
zlabel('n yavg(x,t)')
view([150 25])

function [c,f,s] = N_yavg(x,t,u,dudx)
global l_h Dp_Ds Gp_Ds Gw_Ds

[cSolute, dcSolutedX] = solute_withDerivative(x,t);

VDP = Gp_Ds*dcSolutedX/cSolute; Vs = Gw_Ds*dcSolutedX/cSolute;

calDp = Dp_Ds + (2/105)*l_h^-2*Dp_Ds^-1*Vs*Vs;

c = 1;
f = -VDP*u + calDp*dudx;
s = 0;
end
%----------------------------------------------
function u0 = initial_conditions(x)
global erf_denom
%u0 = erf(x/erf_denom);
u0 = 1;
%u0 = 1 - erf(x/erf_denom);
end
%----------------------------------------------
function [pl,ql,pr,qr] = boundary_conditions(xl,ul,xr,ur,t)
pl = ul-1; 
ql = 0; 
pr = 0;
qr = 1; 
end
%----------------------------------------------

function [C,dCdx] = solute_withDerivative(xx,tt)
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