% solves for y-averaged colloid concentration using 1D approximate equation

global B c0 l_h Dp_Ds 
B = 0.5; c0 = 1.0; l_h = 10; Dp_Ds = 1e-4;

Peclet = 0.334;
PsiW = 0; PsiD = -3;
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
x = linspace(0,1,25);
t = linspace(start_time,2.05,25);

m = 0;
N_yavg_solved = pdepe(m,@N_yavg,@initial_conditions,@boundary_conditions,x,t);

surf(x,t,N_yavg_solved(:,:,1))
xlabel('x')
ylabel('t')
zlabel('n_yavg(x,t)')
view([150 25])

function [c,f,s] = N_yavg(x,t,u,dudx)
global l_h Dp_Ds Gp_Ds Gw_Ds
A1 = (1/105)*(1/(l_h*Dp_Ds)^2)*(7*Gp_Ds*Gw_Ds-4*Gw_Ds^2);
A2 = (1/105)*(1/(l_h*Dp_Ds)^2)*(7*Gp_Ds*Gw_Ds-2*Gw_Ds^2);

[cSolute, dcSolutedX, d2cSolutedX2] = solute_withDerivatives(x,t);
Up_Ds = Gp_Ds*dcSolutedX/cSolute + A1*(dcSolutedX/cSolute^3)*(cSolute*d2cSolutedX2-dcSolutedX^2);
calDp_Ds = Dp_Ds - A2*(dcSolutedX/cSolute)^2;

%Up_Ds = 0;
%calDp_Ds = Dp_Ds;

c = 1/calDp_Ds;
f = dudx;
s = -Up_Ds*dudx/calDp_Ds;
end
%----------------------------------------------
function u0 = initial_conditions(x)
global erf_denom
u0 = erf(x/erf_denom);
end
%----------------------------------------------
function [pl,ql,pr,qr] = boundary_conditions(xl,ul,xr,ur,t)
pl = ul; 
ql = 0; 
pr = 0;
qr = 1; 
end
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
val=C;
C = C + (1-B)*2*(1-cos(lambda))/lambda*sin(lambda*xx)*exp(-lambda^2*tt);
res=sqrt(sum(sum((C-val).^2))./sum(sum(C.^2)));
dCdx = dCdx + (1-B)*2*(1-cos(lambda))*cos(lambda*xx)*exp(-lambda^2*tt);
d2Cdx2 = d2Cdx2 -(1-B)*2*(1-cos(lambda))*lambda*sin(lambda*xx)*exp(-lambda^2*tt);
i=i+1;
end
C=C+B;
end