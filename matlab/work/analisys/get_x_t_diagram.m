function res1 = get_x_t_diagram;

dt = 1e-10;

w = 1.0e9;
T_mod = 2*pi/w;

filter = ones(9,12)/9/12;


% 
% path = 'G:\_results\mobile_16-48e14_10e12_01_Ly_12_01_08\';
% load('G:\_results\mobile_16-48e14_10e12_01_Ly_12_01_08\mobile_16-48e14_10e12_01_Ly_03_Ti.mat', 'param');
% path4save = 'G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\_analisys\';

% path = 'G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\';
% load('G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\ion_motion_16-48_10e12_01Ly_periodic.mat', 'param');



% path = 'G:\_results\ion_motion_16-48e14_2e11_01Ly_periodic_29_12\';
% path4wr = 'G:\_results\ion_motion_16-48e14_2e11_01Ly_periodic_29_12\_movie\';
% load('G:\_results\ion_motion_16-48e14_2e11_01Ly_periodic_29_12\ion_motion_16-48e14_2e11_01Ly_periodic.mat', 'param');

% path = 'I:\_results\ion_motion_16-48_2e11_10Ly_periodic_06_03\';
% load('I:\_results\ion_motion_16-48_2e11_10Ly_periodic_06_03\ion_motion_16-48_2e11_10Ly_periodic.mat', 'param');

% path = 'I:\_results\ion_motion_16-48_6e11_10Ly_periodic_11_03\';
% load('I:\_results\ion_motion_16-48_6e11_10Ly_periodic_11_03\ion_motion_16-48_6e11_10Ly_periodic.mat', 'param');
% 
% path = 'I:\_results\ion_motion_16-48_6e10_10Ly_periodic_14_03\';
% load('I:\_results\ion_motion_16-48_6e10_10Ly_periodic_14_03\ion_motion_16-48_6e10_10Ly_periodic.mat', 'param');

% path = 'G:\_results\moveless_16-48e14_thin_beam10e12_07_21\';
% load('G:\_results\moveless_16-48e14_thin_beam10e12_07_21\moveless_16-48e14_thin_beam10e12_07_21.mat', 'param');

% path = 'G:\_results\moveless_48-16e14_10e12_01Ly_07_27\';
% load('G:\_results\moveless_48-16e14_10e12_01Ly_07_27\moveless_48-16e14_thin_beam10e12_07_27.mat', 'param');

path = 'G:\_results\ion_motion_16-48_10e12_01Ly_periodic_02_16_detailed\';
load('G:\_results\ion_motion_16-48_10e12_01Ly_periodic_02_16_detailed\ion_motion_16-48_10e12_01Ly_periodic.mat', 'param');


geometry = param.geometry;
bc = param.bc;
geometry = calc_grid_step(geometry,bc);
rho_back_struct = param.rho_back_struct;

    rho_back = load_rho_back(rho_back_struct, geometry, bc)*1.6e-19;


var = 'rho_sp';
var1 = 'rho';

% e_abs_full = zeros(2048,1);
% rho_ions_full = zeros(2048,1);

% rho_ions_x_t = zeros(256,1:32*60);
grid_position = round(0.27/60*2048);

% e_abs_x_920 = zeros(256,1);

tstart = 0e-7;
tend = 3.9e-7;
dt = 1e-9/5;

t = tstart:dt:tend;

% rho_ions_full = zeros(2048,length(t));

% e_abs_full = zeros(2048,length(t));
% rho_electrons_full = zeros(2048,length(t));
ex_full = zeros(2048,length(t));
i = 1;

for t = tstart:dt:tend
    t
    
%     load([path,'_analisys\',num2str(i,'%3.3d'),'_1.dat'],'e_abs_t','-mat');
%     e_abs_full = cat(2,e_abs_full, squeeze(e_abs_t(128,:,:)));
    
%     load([path,'_analisys\',num2str(i,'%3.3d'),'_1_ions.dat'],'rho_ions','-mat');
%     rho_ions = rho_ions-rho_back;

        try,
            
        load('-mat',strcat(path,'rho_electrons_',num2str(t),'.dat'),'rho_sp') 
        rho_el = rho_sp;

        load('-mat',strcat(path,'rho_ions_',num2str(t),'.dat'),var) 
        rho_ions = rho_sp;     
        load('-mat',strcat(path,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;
        catch,
            
            ttt = t+dt/20;
                load('-mat',strcat(path,'rho_electrons_',num2str(ttt),'.dat'),'rho_sp') 
        rho_el = rho_sp;

        load('-mat',strcat(path,'rho_ions_',num2str(ttt),'.dat'),var) 
        rho_ions = rho_sp;     
        load('-mat',strcat(path,'rho_light_electrons_',num2str(ttt),'.dat'),var) 
        rho_beam = rho_sp;  

        end
        
        %         rho = rho_sp;
        rho = rho_el+rho_ions+rho_beam;
%         rho = rho_el+rho_back+rho_beam;
% % 
% %         

        fi = field_3(rho,geometry,bc);
% 
        [ex ey] = e_from_fi(fi,geometry,bc);
%         e_abs = (ex.^2 + ey.^2);
        
        
%         rho_ions = imfilter(rho_ions,filter) - rho_back;    
%         rho_ions = rho_sp - rho_back; 
        
%         e_abs_full(:,i) = e_abs(128,:);
        ex_full(:,i) = ex(128,:);
%         rho_ions_full(:,i) = rho_ions(128,:);
% rho_el = imfilter(rho_el,filter) + rho_back;
%     rho_electrons_full(:,i) = rho_el(128,:);
        i = i+1;
 

    
        
    
%     load([path,'_analisys\',num2str(i,'%3.3d'),'_1.dat'],'e_abs_t','-mat');
%     e_abs_x_920 = cat(2,e_abs_x_920, squeeze(e_abs_t(:,920,:)));
    
    
end

% res1 = e_abs_x_920(:,2:end);
% res1 = rho_ions_full;
res1 = ex_full(:,1:end);
% res1 = rho_electrons_full;
    