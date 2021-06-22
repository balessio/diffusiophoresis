% solves for y-averaged colloid concentration using 1D approximate equation

global B c0 Dp_Ds 
B = 0.1; c0 = 1.0; Dp_Ds = 1e-1; h_l = 0.05; w_l = 0.05;

%w_l = 0.0

Peclet = 0.334;
%global PsiW
PsiW = -4;
PsiD = -3;
D_pos_D_Na = 1.0; D_neg_D_Na = 1.526;
global Gw_Ds Gp_Ds
Gw_Ds = Peclet*(((D_pos_D_Na-D_neg_D_Na)/(D_pos_D_Na+D_neg_D_Na))*PsiW+PsiW^2/8);
Gp_Ds = Peclet*(((D_pos_D_Na-D_neg_D_Na)/(D_pos_D_Na+D_neg_D_Na))*PsiD+PsiD^2/8);
global modDiffConst
modDiffConst = braket_ff(h_l,w_l)
%----------------------------------------------

start_time = 0.001;
global erf_denom
erf_denom = sqrt(4*start_time);
if start_time==0.0
   erf_denom = 0.01; 
end

%x = [linspace(0,0.09,100) linspace(0.1,1,10)];
x = linspace(0,1,100);
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

Vs = (x>0.06)*Gw_Ds*dudx(2)/u(2);
%Vs = Gw_Ds*dudx(2)/u(2);

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
u0 = [exp(-(x-0.06)^2/0.002);B*c0 + c0*(1-B)*erf(x/erf_denom)];
end
%----------------------------------------------
function [pl,ql,pr,qr] = boundary_conditions(xl,ul,xr,ur,t)
global B c0
%pl = [0;ul(2)-B*c0];
pl = [ul(1);ul(2)-B*c0];%[ul(1)-1;ul(2)-B*c0]; 
%ql = [1;0];
ql = [0;0]; 
pr = [0;0];
qr = [1;1]; 
end
%----------------------------------------------
%----------------------------------------------
function lambda_k = lambda_coefficient(k)
	lambda_k = (2*k+1)*pi;
end
%----------------------------------------------
function a_k = a_coefficient(k,lambda_k,w_h)
	a_k = 8*(-1)^k*lambda_k^-3/cosh(0.5*lambda_k*w_h);
end
%----------------------------------------------
function b_mk = b_coefficient(m,k,alpha_const,lambda_m,lambda_k,a_m,b0_mk,h_l,w_l)
	w_h = w_l/h_l;
	b_mk = 4*(alpha_const*(-1)^(k+1)*a_m*lambda_k...
		*cosh(0.5*lambda_m*w_h)/((lambda_m*w_h)^2+lambda_k^2) ...
		+ (1-alpha_const/6)*b0_mk)...
	/((lambda_m/h_l)^2+(lambda_k/w_l)^2);
end
%----------------------------------------------
function b0_mk = b0_coefficient(m,k,lambda_m,lambda_k)
	b0_mk = 4*(-1)^(m+k)/lambda_m/lambda_k;
end
%----------------------------------------------
% we calculate the dispersion in four pieces: K1, K2, K3, K4
function diffusivity_constant = braket_ff(h_l,w_l)
sumToIndex = 10;
w_h = w_l/h_l;
k = 0; alpha_const = 0;
while k <= sumToIndex
	lambda_k = lambda_coefficient(k);
	alpha_const = alpha_const + (2/lambda_k)^5*tanh(0.5*lambda_k*w_h);
	k = k+1;
end
alpha_const = 6/(1-6*alpha_const/w_h);
fprintf('alpha_const=')
disp(alpha_const)

K1 = alpha_const^2*h_l^2/7560;
fprintf('K1=')
disp(K1)

k = 0; K2 = 0;
while k <= sumToIndex
	lambda_k = lambda_coefficient(k); a_k = a_coefficient(k,lambda_k,w_h);
	K2 = K2 + (-1)^(k+1)*a_k*sinh(0.5*lambda_k*w_h)...
	*(lambda_k^-2 + 60*lambda_k^-4 - 720*lambda_k^-6);
	k = k+1;
end
K2 = alpha_const^2*h_l^2/w_h*K2/90;
fprintf('K2=')
disp(K2)

m = 0; K3 = 0;
while m <= sumToIndex
	lambda_m = lambda_coefficient(m); a_m = a_coefficient(m,lambda_m,w_h);
	k = 0;
	K3sum = 0;
	while k <= sumToIndex
		lambda_k = lambda_coefficient(k); b0_mk = b0_coefficient(m,k,lambda_m,lambda_k);
		b_mk = b_coefficient(m,k,alpha_const,lambda_m,lambda_k,a_m,b0_mk,h_l,w_l);
		K3sum = K3sum + (-1)^(m+k)*b_mk*(lambda_m^2-12)/lambda_m^3/lambda_k;
		k = k+1;
	end
	K3 = K3 + K3sum;
	m = m+1;
end
K3 = 2*alpha_const*K3/3;
fprintf('K3=')
disp(K3)

m = 0; K4 = 0;
while m <= sumToIndex
	lambda_m = lambda_coefficient(m); a_m = a_coefficient(m,lambda_m,w_h);
	k=0; K4sum = 0;
	while k <= sumToIndex
		lambda_k = lambda_coefficient(k); b0_mk = b0_coefficient(m,k,lambda_m,lambda_k);
		b_mk = b_coefficient(m,k,alpha_const,lambda_m,lambda_k,a_m,b0_mk,h_l,w_l);
		K4sum = K4sum + a_m*b_mk*((-1)^(k+1)*lambda_k*(lambda_k^2/w_h+lambda_m^2*w_h)^(-1)...
		*cosh(0.5*lambda_m*w_h) - 4*(-1)^(m+1)*b0_mk*sinh(0.5*lambda_m*w_h)/lambda_m^2);
		k = k+1;
	end
	K4 = K4 + K4sum;
	m = m+1;
end
K4 = 2*alpha_const*K4/w_l^2;
fprintf('K4=')
disp(K4)

diffusivity_constant = K1+K2+K3+K4;
end
