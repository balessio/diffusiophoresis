% email balessio@princeton.edu for correspondence

% solves for solute concentration for a given y-position
% we just do one solute here (c)

% outputs two cells 
% x, y, t

% {c_t1(x,y)   ...  c_t1_soluteN(x,y)
% c_t2_solute1(x,y)  ...  c_t2_soluteN(x,y)
% c_t3_solute1(x,y)  ...  c_t3_soluteN(x,y)
% ...
% c_tfinal_solute1(x,y)}   ...  c_tfinal_soluteN(x,y)

% inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global parameter_file_upload
parameter_file_unload = matfile(parameter_file_upload);

global Peclet L_h delta_x delta_y delta_t start_time duration
global z_i D_i PsiWref cWref

Peclet = parameter_file_unload.Peclet;
L_h = parameter_file_unload.L_h;
delta_x = parameter_file_unload.delta_x;
delta_y = parameter_file_unload.delta_y;
delta_t = parameter_file_unload.delta_t;
start_time = parameter_file_unload.start_time;
duration = parameter_file_unload.duration;
z_i = parameter_file_unload.z_i;
D_i = parameter_file_unload.D_i;
PsiWref = parameter_file_unload.PsiWref;
cWref = parameter_file_unload.cWref;

output_folder = parameter_file_unload.output_folder;

global  numY 
numY = 1/(delta_y*L_h)-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

x = 0:delta_x:1;
t = start_time:delta_t:duration;

m = 0;
options=odeset('NonNegative',[]);

tic;
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t, options);
toc;

y = [flip(-(delta_y:delta_y:1/L_h)) 0:delta_y:1/L_h];

PsiWref = ones(length(x),1)*PsiWref;

XYTsolute_cell=cell(1,3); XYTsolute_cell{1,1}=x;
XYTsolute_cell{1,2}=y; XYTsolute_cell{1,3}=t;
solute_cell=cell(length(t),length(z_i));
middle_Y_index = int8((length(y)+1)/2);
for i = 1:length(t)
    for ion=1:length(z_i)
        stack = int8((ion-1)*numY); % each ion is stacked in the solution vector
        c_frame = zeros(length(x),length(y));
        c_frame(:,1) = sol(i,:,numY+stack);
        c_frame(:,length(y)) = sol(i,:,numY+stack);
        c_frame(:,middle_Y_index) = sol(i,:,1+stack);
        for j = 1:middle_Y_index-2
            c_frame(:,middle_Y_index-j) = sol(i,:,j+stack);
            c_frame(:,middle_Y_index+j) = sol(i,:,j+stack);
        end
        solute_cell{i,ion} = c_frame;
    end
end

save(join([output_folder,'XYTsolute_cell.mat'],''), 'XYTsolute_cell');
save(join([output_folder,'solute_cell.mat'],''), 'solute_cell');

figure(1);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
ax(1) = subplot(2,3,1);
surf(x,t,sol(:,:,1))
xlabel('Distance X')
ylabel('Time T')
zlabel('solute concentration c1(X,Y=0,T)')
ax(2) = subplot(2,3,2);
surf(x,t,sol(:,:,1+numY))
xlabel('Distance X')
ylabel('Time T')
zlabel('solute concentration c2(X,Y=0,T)')
ax(3) = subplot(2,3,3);
surf(x,t,sol(:,:,1+2*numY))
xlabel('Distance X')
ylabel('Time T')
zlabel('solute concentration c3(X,Y=0,T)')
ax(1) = subplot(2,3,4);
surf(x,t,sol(:,:,numY-1))
xlabel('Distance X')
ylabel('Time T')
zlabel('solute concentration c1(X,Y=h/L,T)')
ax(2) = subplot(2,3,5);
surf(x,t,sol(:,:,numY-1+numY))
xlabel('Distance X')
ylabel('Time T')
zlabel('solute concentration c2(X,Y=h/L,T)')
ax(3) = subplot(2,3,6);
surf(x,t,sol(:,:,numY-1+2*numY))
xlabel('Distance X')
ylabel('Time T')
zlabel('solute concentration c3(X,Y=h/L,T)')

saveas(gcf,join([output_folder,'3D_plot.fig'],''));

%----------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,dudx) %#ok<*INUSL> % Equation to solve
global numY delta_y Peclet PsiWref cWref D_i z_i
Vs_num1=0; Vs_denom1=0; Vs_num2=0; Vs_denom2=0;
for ion=1:length(z_i)
   Vs_num1 = Vs_num1 + D_i(ion)*z_i(ion)*dudx(ion*numY);
   Vs_denom1 = Vs_denom1 + D_i(ion)*z_i(ion)^2*u(ion*numY);
   Vs_num2 = Vs_num2 + z_i(ion)^2*dudx(ion*numY);
   Vs_denom2 = Vs_denom2 + z_i(ion)^2*u(ion*numY);
end
% solve_PsiW = @(PsiW) F_minus_Fref_W(PsiW,u(numY),...
%     u(2*numY),u(3*numY));
% PsiW_sol = fsolve(solve_PsiW,PsiWref);

c_array = zeros(length(x),3);
c_array(:,1) = u(numY);
c_array(:,2) = u(2*numY);
c_array(:,3) = u(3*numY);
    
solve_PsiW = @(PsiW) F_minus_Fref(PsiW,c_array,PsiWref,cWref,z_i);
PsiW_sol = fsolve(solve_PsiW,PsiWref);

Vs = Peclet*(Vs_num1*PsiW_sol/Vs_denom1 + Vs_num2*PsiW_sol^2/Vs_denom2/8);
c = ones(length(z_i)*numY,1);
f_array_make = zeros(length(z_i)*numY,1);
s_array_make = zeros(length(z_i)*numY,1);
for i = 1:numY
    Vfx = 0.5*Vs*(3*(i*delta_y)^2-1);
    i_plus_1 = i+1; i_minus_1 = i-1;
    if i==1
        i_minus_1 = i;
    end
    if i==numY
        i_plus_1 = i;
    end
    x_em_num=0; x_em_denom=0;
    % six different summation terms from the electromigration expression in
    % the y-direction
    y_em_s1=0; y_em_s2=0; y_em_s3=0; y_em_s4=0; y_em_s5=0; y_em_s6=0;
    for ion=1:length(z_i)
        stack = int8((ion-1)*numY);
        DCDY = (u(i_plus_1+stack)-u(i_minus_1+stack))/(2*delta_y);
        D2CDY2 = (u(i_plus_1+stack)-2*u(i+stack)+u(i_minus_1+stack))/delta_y^2;
        x_em_num = x_em_num + z_i(ion)*D_i(ion)*dudx(i+stack);
        x_em_denom = x_em_denom + z_i(ion)^2*D_i(ion)*u(i+stack);
        y_em_s1 = y_em_s1 + z_i(ion)*D_i(ion)*D2CDY2;
        y_em_s2 = y_em_s2 + z_i(ion)*D_i(ion)*DCDY;
        y_em_s3 = y_em_s3 + z_i(ion)^2*D_i(ion)*u(i+stack);
        y_em_s5 = y_em_s5 + z_i(ion)^2*D_i(ion)*DCDY;
    end
    y_em_s4 = y_em_s2;
    y_em_s6 = y_em_s3;
    x_em = x_em_num/x_em_denom;
    for ion = 1:length(z_i)
        stack = int8((ion-1)*numY); % each ion is stacked in the solution vector
        DCDY = (u(i_plus_1+stack)-u(i_minus_1+stack))/(2*delta_y);
        D2CDY2 = (u(i_plus_1+stack)-2*u(i+stack)+u(i_minus_1+stack))/delta_y^2;
        y_partFlux_c = -Vs*(alpha(i_plus_1)*u(i_plus_1+stack)...
            - alpha(i_minus_1)*u(i_minus_1+stack));
        f_array_make(i+stack) = D_i(ion)*dudx(i+stack)...
            - z_i(ion)*D_i(ion)*u(i+stack)*x_em...
            - Vfx*u(i+stack) + y_partFlux_c;
        y_partSource_c = Vs*(alpha(i_plus_1)*dudx(i_plus_1+stack)...
            - alpha(i_minus_1)*dudx(i_minus_1+stack))...
            - z_i(ion)*D_i(ion)... 
            *((u(i+stack)*y_em_s1+DCDY*y_em_s2)*y_em_s3...
            - u(i+stack)*y_em_s4*y_em_s5)/y_em_s6^2;
        s_array_make(i+stack) = D_i(ion)*D2CDY2 + y_partSource_c;
    end
end
f = f_array_make;
s = s_array_make;
end
%----------------------------------------------
function u0 = pdex1ic(x) % Initial conditions
global numY start_time
u0 = zeros(3*numY,1);
erf_denom = sqrt(4*start_time);
if start_time==0.0
   erf_denom = 0.01;
end
for i = 1:numY
    u0(i) = 1.0*erf(x/erf_denom);
    u0(i+numY) = 1.0;
    u0(i+2*numY) = 1.0*(1.0-erf(x/erf_denom));
end
end
%----------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t) % Boundary conditions
global numY
pl = zeros(3*numY,1);
for i = 1:numY
   pl(i) = ul(i) - 0.0;
   pl(i+numY) = ul(i+numY) - 1.0;
   pl(i+2*numY) = ul(i+2*numY) - 1.0;
end
ql = zeros(3*numY,1);
pr = zeros(3*numY,1);
qr = ones(3*numY,1); % dudx = 0, this term just needs to be nonzero
end
%----------------------------------------------
function solve_for_Psi = F_minus_Fref(Psi_,c_,Psi_ref,c_ref,z_i) % solve PsiW with refs
s = size(c_);
c_ref_array = zeros(s);
for ion = 1:length(z_i)
   c_ref_array(:,ion) = ones(s(1),1)*c_ref(ion);
end
solve_for_Psi = FF(Psi_,c_,z_i)-FF(Psi_ref,c_ref_array,z_i);
end
function Fpsifunc = FF(Psi_,c_,z_i) % solve for Psi with references
s = size(c_);
Fpsifunc = zeros(s(1),1);
for i=1:length(z_i)
    Fpsifunc = Fpsifunc + c_(:,i).*(exp(-z_i(i).*Psi_) - 1);
end
Fpsifunc = 2*Fpsifunc;
end
%----------------------------------------------
function alpha_return = alpha(j) % coefficient dependent on Y
global delta_y L_h
alpha_return = -0.25*j*((j*delta_y*L_h)^2-1);
end
%----------------------------------------------
