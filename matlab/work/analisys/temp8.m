function [res1 res2 res3] = temp8

path = 'G:\_results\mobile_20-44e14_05e11_01Ly_04.05\';
time = 6.6e-7;
specie = 'rho_ions_';

load('G:\_results\mobile_20-44e14_05e11_01Ly_04.05\mobile_20-44e14_05e11_01Ly.mat', 'param');
geometry = param.geometry;
bc = param.bc;
rho_back_struct = param.rho_back_struct;

geometry = calc_grid_step(geometry, bc);
const.eq = 1.6e-19;

rho_back = load_rho_back(rho_back_struct, geometry, bc)*const.eq;

h = ones(10)/100;

t_start = 3.5e-7;
t_end = 3.74e-7;
dt = 5e-10;

N = length(t_start:dt:t_end);

e_abs = zeros(256,2048);
ex_av = zeros(256,2048);
ey_av = zeros(256,2048);

for time = t_start:dt:t_end
% 
time
% time = 6.39e-7;

load('-mat', strcat(path, 'rho_ions_', num2str(time),'.dat'), 'rho_sp');

rho_ions = rho_sp;

load('-mat', strcat(path, 'rho_electrons_', num2str(time),'.dat'), 'rho_sp');

rho_el = rho_sp;
load('-mat', strcat(path, 'rho_light_electrons_', num2str(time),'.dat'), 'rho_sp');

rho_beam = rho_sp;

rho = (rho_el + rho_ions + rho_beam);

fi = field_3(rho, geometry, bc);

[ex ey] = e_from_fi(fi, geometry, bc);

e_abs = (ex.^2 + ey.^2).^0.5 + e_abs;
% ex_av = ex_av + ex;
% ey_av = ey_av + ey;
% 
end
% 
e_abs = e_abs/N;
% ex_av = ex_av/N;
% ey_av = ey_av/N;

rho = rho_ions-rho_back;

rho = rho/const.eq;

rho = imfilter(rho,h);


f1 = figure('Position', [10 10 600 300], 'Color', 'w');

% figure, imshow(ex_av,[])
% figure, imshow(ey_av,[])
% figure, imshow(e_abs,[])
% figure, imshow((ex_av.^2 + ey_av.^2).^0.5,[])



rho_beam = rho_beam/const.eq;

clim1 = [-3.0e14 4.0e14];
clim2 = [0 8e3];
clim3 = [-3e12 0];

max(max(e_abs))

a1 = axes('Parent', f1);
set(a1, 'Unit', 'normalized', 'Position', [0.15 0.15 0.8 0.8]);
im1 = image(e_abs(:,683:683*2-1),'Parent',a1, 'CDataMapping', 'scaled');
set(a1, 'xTick', [10 340 680], 'xTickLabel', ['0.2'; '0.3'; '0.4'],...
                'yTick', [10 128 256], 'yTickLabel', [' 0.0'; '0.05'; ' 0.1'], 'YDir', 'normal','Clim',clim2);
%              'Clim', clim1); 

colormap('gray')
% figure, surf(rho(:,700:1024+100))

res1 = ex_av;
res2 = ey_av;
res3 = e_abs;
