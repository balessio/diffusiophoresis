% solves for y-averaged colloid concentration using 1D approximate equation

global B c0 h_l Dp_Ds 
B = 0.1; c0 = 1.0; h_l = 1/20; Dp_Ds = 1e-3; % h = total width

Peclet = 0.334;
%global Psi
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
x = linspace(0,1,200);
t = linspace(start_time,2.05,205);

m = 0;
N_yavg_solved = pdepe(m,@N_yavg,@initial_conditions,@boundary_conditions,x,t);

surf(x,t,N_yavg_solved(:,:,1))
xlabel('x')
ylabel('t')
zlabel('n yavg(x,t)')
view([150 25])


function [c,f,s] = N_yavg(x,t,u,dudx)
global h_l Dp_Ds Gp_Ds Gw_Ds

VDP = Gp_Ds*dudx(2)/u(2);

if x>0
	Vs = Gw_Ds*dudx(2)/u(2);
else
	Vs = 0;
end 


modDiff2D = h_l^2/210;

calDp = Dp_Ds + modDiff2D*Vs^2/Dp_Ds;

calDs = 1 + modDiff2D*Vs^2;

c = [1;1];
f = [-VDP*u(1) + calDp*dudx(1); calDs*dudx(2)];
s = [0;0];
end
%----------------------------------------------
function u0 = initial_conditions(x)
global erf_denom B c0
%u0 = [1;B*c0 + c0*(1-B)*erf(x/erf_denom)];
%u0 = [exp(-(x-0.2)^2/0.01);B*c0 + c0*(1-B)*erf(x/erf_denom)];
u0 = [exp(-(x-0.5)^2/0.01);B*c0 + c0*(1-B)*erf(x/erf_denom)];
%u0 = [exp(-(x-0.06)^2/0.002);B*c0 + c0*(1-B)*erf(x/erf_denom)];
end
%----------------------------------------------
function [pl,ql,pr,qr] = boundary_conditions(xl,ul,xr,ur,t)
global B c0 Gw_Ds Gp_Ds

%dcdx = 0; k = 0;
%while k<15
%	dcdx = dcdx + 2*(1-cos(pi*(2*k+1)/2))*exp(-(pi*(2*k+1)/2)^2*t);
%	k = k+1;
%end
%dcdx = (1-B)*dcdx;
%vp = (Gw_Ds+Gp_Ds)*log(dcdx);

%pl = [0;ul(2)-B*c0]; ql = [1;0]; % no flux
pl = [ul(1);ul(2)-B*c0]; ql = [0;0]; % n=0
%pl = [ul(1)*vp;ul(2)-B*c0]; ql = [1;0]; 

pr = [0;0]; qr = [1;1]; 
end
%----------------------------------------------
