% email balessio@princeton.edu for correspondence

% solves for solute concentration for a given y-position
% we just do one solute here (c)

% outputs a cell of solute concentration matrices and vectors x, y, t
% {c_t1(x,y), x, y, t
% c_t2(x,y)
% c_t3(x,y)
% ...
% c_tfinal(x,y)}

% inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global parameter_file_upload
parameter_file_unload = matfile(parameter_file_upload);

global Peclet B c0 L_h start_time duration delta_x delta_y delta_t

Peclet = parameter_file_unload.Peclet;
B = parameter_file_unload.B;
c0 = parameter_file_unload.c0;
L_h = parameter_file_unload.L_h;
PsiW = parameter_file_unload.PsiW;
D_pos_D_Na = parameter_file_unload.D_pos_D_Na;
D_neg_D_Na = parameter_file_unload.D_neg_D_Na;
start_time = parameter_file_unload.start_time;
duration = parameter_file_unload.duration;
delta_x = parameter_file_unload.delta_x;
delta_y = parameter_file_unload.delta_y;
delta_t = parameter_file_unload.delta_t;

output_folder = parameter_file_unload.output_folder;

global Gw_Ds
Gw_Ds = Peclet*(((D_pos_D_Na-D_neg_D_Na)/(D_pos_D_Na+D_neg_D_Na))*PsiW+(PsiW^2)/8);
global numY
numY = 1/(delta_y*L_h) - 1; % spans y over half width bc symmetry, excludes edges

if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

save(join([output_folder,'log.mat'],''),...
'delta_x', 'delta_t',...
'start_time', 'duration', 'B', 'c0', 'L_h', 'Gw_Ds');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the grid
% start at t=0.001, which corresponds to pure diffusion of form
% erf(x/sqrt(4*0.001)) (see Ault et al equation 5)
x = 0:delta_x:1;

t = start_time:delta_t:duration;
m = 0;
options=odeset('NonNegative',[]);

tic;
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t, options);
toc;

y = [flip(-(delta_y:delta_y:1/L_h)) 0:delta_y:1/L_h];

XYTsolute_cell=cell(1,3); XYTsolute_cell{1,1}=x;
XYTsolute_cell{1,2}=y; XYTsolute_cell{1,3}=t;
solute_cell=cell(length(t),1);
middle_Y_index = int8((length(y)+1)/2);
for i = 1:length(t)
    c_frame = zeros(length(x),length(y));
    c_frame(:,1) = sol(i,:,numY);
    c_frame(:,length(y)) = sol(i,:,numY);
    c_frame(:,middle_Y_index) = sol(i,:,1);
    for j = 1:middle_Y_index-2
        c_frame(:,middle_Y_index-j) = sol(i,:,j);
        c_frame(:,middle_Y_index+j) = sol(i,:,j);
    end
    solute_cell{i,1} = c_frame;
end

save(join([output_folder,'XYTsolute_cell.mat'],''), 'XYTsolute_cell');
save(join([output_folder,'solute_cell.mat'],''), 'solute_cell');

figure(1);
surf(x,t,sol(:,:,1))
xlabel('Distance X')
ylabel('Time T')
zlabel('colloid concentration c(X,Y=0,T)')
saveas(gcf,join([output_folder,'3D_plot.png'],''));

%----------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,dudx) %#ok<*INUSL> % Equation to solve
global L_h Gw_Ds numY delta_y
Vs = -Gw_Ds*dudx(numY)/u(numY);
c = zeros(numY,1) + 1;
f_array_make = zeros(numY,1);
s_array_make = zeros(numY,1);
for i = 1:numY
    i_plus_1 = i+1; i_minus_1 = i-1;
    y_m = i_minus_1*delta_y; y_p = i_plus_1*delta_y;
    if i==1
       y = 0; 
    elseif i==numY
       y = 1/L_h; 
    else
        y = i*delta_y;
    end
    Vfx = 0.5*Vs*(3*(y*L_h)^2-1);
    
    if i==1
        i_minus_1 = i; y_m = 0;
    elseif i==numY
        i_plus_1 = i; y_p = 1/L_h;
    end

    y_partFlux_c = -(0.5/delta_y)*0.5*Gw_Ds...
        *((alpha(y_p)*u(i_plus_1)...
        -alpha(y_m)*u(i_minus_1))*dudx(numY)/u(numY));
	f_array_make(i) = dudx(i)-Vfx*u(i) + y_partFlux_c;
    y_partSource_c = -(0.5/delta_y)*0.5*Gw_Ds*y*((y*L_h)^2-1)...
        *(-(alpha(y_p)*dudx(i_plus_1)+alpha(y_m)*dudx(i_minus_1))...
        *dudx(numY)/u(numY)...
        +(alpha(y_p)*u(i_plus_1)-alpha(y_m)*u(i_minus_1))...
        *dudx(numY)^2/u(numY)^2);
    s_array_make(i) = (1/delta_y^2)...
        *(u(i_plus_1)-2*u(i)+u(i_minus_1)) + y_partSource_c;
end
f = f_array_make; % might need to transpose by: f'
s = s_array_make;
end
%----------------------------------------------
function u0 = pdex1ic(x) % Initial conditions
global c0 B numY start_time
u0 = zeros(numY,1);
erf_denom = sqrt(4*start_time);
if start_time==0.0
   erf_denom = 0.01; 
end
for i = 1:numY
    u0(i) = B*c0 + c0*(1-B)*erf(x/erf_denom);
end
end
%----------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t) % Boundary conditions
global B c0 numY 
pl = zeros(numY,1);
for i = 1:numY
   pl(i) = ul(i)-B*c0;
end
ql = zeros(numY,1);
pr = zeros(numY,1);
qr = zeros(numY,1) + 1; % dudx = 0, this term just needs to be nonzero
end
%----------------------------------------------
function alpha_coeff = alpha(yy)
global L_h
alpha_coeff = yy*((yy*L_h)^2-1);
end
%----------------------------------------------