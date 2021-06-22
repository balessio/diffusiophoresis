% solves for y-averaged colloid concentration using 1D approximate equation

global B c0 Dp_Ds 
%global W
B = 0.1; c0 = 1.0; Dp_Ds = 1e-2; H = 0.05; W = 0.06;

Peclet = 0.334;
%global PsiW
PsiW = -4;
PsiD = -3;
D_pos_D_Na = 1.0; D_neg_D_Na = 1.526;
global Gw_Ds Gp_Ds
Gw_Ds = Peclet*(((D_pos_D_Na-D_neg_D_Na)/(D_pos_D_Na+D_neg_D_Na))*PsiW+PsiW^2/8);
Gp_Ds = Peclet*(((D_pos_D_Na-D_neg_D_Na)/(D_pos_D_Na+D_neg_D_Na))*PsiD+PsiD^2/8);
global modDiffConst
modDiffConst = braket_ff(H,W)
%----------------------------------------------

start_time = 0.001;
global erf_denom
erf_denom = sqrt(4*start_time);
if start_time==0.0
   erf_denom = 0.01; 
end

%x = [linspace(0,0.09,100) linspace(0.1,1,10)];
x = linspace(0,1,200);
t = linspace(start_time,1.05,50);

m = 0;
N_yavg_solved = pdepe(m,@N_yavg,@initial_conditions,@boundary_conditions,x,t);

surf(x,t,N_yavg_solved(:,:,1))
xlabel('x')
ylabel('t')
zlabel('n yavg(x,t)')
view([150 25])

function [c,f,s] = N_yavg(x,t,u,dudx)
global l_h Dp_Ds Gp_Ds Gw_Ds modDiffConst

VDP = Gp_Ds*dudx(2)/u(2);

%Vs = (x>0.06)*Gw_Ds*dudx(2)/u(2);
Vs = Gw_Ds*dudx(2)/u(2);

calDp = Dp_Ds + modDiffConst*Dp_Ds^-1*Vs^2;
calDs = 1 + modDiffConst*Vs^2;

c = [1;1];
f = [-VDP*u(1) + calDp*dudx(1); calDs*dudx(2)];
s = [0;0];
end
%----------------------------------------------
function u0 = initial_conditions(x)
global erf_denom B c0
%u0 = [1;B*c0 + c0*(1-B)*erf(x/erf_denom)];
%u0 = [exp(-(x-0.06)^2/0.002);B*c0 + c0*(1-B)*erf(x/erf_denom)];
u0 = [exp(-(x-0.2)^2/0.01);B*c0 + c0*(1-B)*erf(x/erf_denom)];
end
%----------------------------------------------
function [pl,ql,pr,qr] = boundary_conditions(xl,ul,xr,ur,t)
global B c0
pl = [0;ul(2)-B*c0]; ql = [1;0];
%pl = [ul(1);ul(2)-B*c0]; ql = [0;0];
pr = [0;0]; qr = [1;1]; 
end
%----------------------------------------------
%----------------------------------------------
function A = A_const(W_H)
	sum_to = 50;
	A = 0; k = 0;
	while k < sum_to
		lambda_k = (k+0.5)*pi; A = A + tanh(lambda_k*W_H)/lambda_k^5;
		k = k+1;
	end
	A = 1.5/(1-6*A/W_H)
end
%----------------------------------------------
function alpha_mp = alpha_coeff(m,p,W_H)
	sum_to = 50;
	alpha_mp = 0; k = 0;
	while k < sum_to
		alpha_mp = alpha_mp + tanh((k+0.5)*pi*W_H)...
		/((k+0.5)*((k+0.5)^2-m^2)*((k+0.5)^2+p^2/W_H^2));
		k = k+1;
	end
	alpha_mp = -16/pi^5/W_H/(-1)^(m+p+1)*alpha_mp/(1+1*(m==0)+1*(p==0));
end
%----------------------------------------------
% we calculate the dispersion in three pieces: K1, K2, K3
function diffusivity_constant = braket_ff(H,W)
sum_to = 10;
W_H = W/H;

A = A_const(W_H);
fprintf('A = ')
disp(A)

K1 = A^2*H^2*2/945;
fprintf('K1 = ')
disp(K1)

K2 = 0; k = 1;
while k < sum_to
	alpha_k0 = alpha_coeff(k,0,W_H);
	K2 = K2 + (-1)^(k+1)*alpha_k0/k^4;
	k = k+1;
end
K2 = A^2*H^2*K2/pi^4;
fprintf('K2 = ')
disp(K2)

K3 = 0; m = 0;
while m < sum_to
	pTerm = 0; p = 0;
	while p < sum_to
		alpha_mp = alpha_coeff(m,p,W_H);
		if (m~=0||p~=0)
			pTerm = pTerm + (1+1*(m==0)+1*(p==0))*alpha_mp^2/(m^2+p^2/W_H^2);
		end
		p = p+1;
	end
	K3 = K3 + pTerm;
	m = m+1;
end
K3 = A^2*H^2*K3/16/pi^2;
fprintf('K3 = ')
disp(K3)

diffusivity_constant = K1+K2+K3;
end
