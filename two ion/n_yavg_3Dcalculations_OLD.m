% solves for y-averaged colloid concentration using 1D approximate equation

global B c0 Dp_Ds 
B = 0.1; c0 = 1.0; Dp_Ds = 1e-1; h_l = 0.05; w_l = 0.1;

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
