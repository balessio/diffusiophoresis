% solves for y-averaged colloid concentration using 1D approximate equation

H = 0.025; %W = 0.1;
global W 

modDiffConst = braket_ff(H,W)

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
