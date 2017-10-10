function res1 = calc_alt_field(amb_field_x, amb_field_y);

dt = 1e-10;

w = 1.0e9;
T_mod = 2*pi/w;

filter = ones(9,12)/9/12;


% 
% path = 'G:\_results\mobile_16-48e14_10e12_01_Ly_12_01_08\';
% load('G:\_results\mobile_16-48e14_10e12_01_Ly_12_01_08\mobile_16-48e14_10e12_01_Ly_03_Ti.mat', 'param');
% path4save = 'G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\_analisys\';

path = 'G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\';
load('G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\ion_motion_16-48_10e12_01Ly_periodic.mat', 'param');


geometry = param.geometry;
bc = param.bc;
geometry = calc_grid_step(geometry,bc);
rho_back_struct = param.rho_back_struct;

var = 'rho_sp';
var1 = 'rho';


% if rho_back_struct.enabled
    rho_back = load_rho_back(rho_back_struct, geometry, bc)*1.6e-19;
% end

    j = 1;
    dt = 1e-9;
%     dt = 5e-10;
    
    t_start = 2.6e-7 ;
    t_end = t_start + 4*T_mod - dt;

    
    e_abs_av = zeros(geometry.ngy, geometry.ngx);
    
    
    
    
    N = length(t_start:dt:t_end)

    
    
            
            for t = t_start:dt:t_end
                
%                 (t_start + k*dt*n_slices):dt:(t_start + (k + 1)*dt*n_slices - dt)
                t
                
               
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
                    
                rho = rho_el + rho_ions + rho_beam;    
                
                fi = field_3(rho,geometry,bc);
                [ex ey] = e_from_fi(fi,geometry,bc);
%                 e_abs = (ex.^2 + ey.^2);
                
%                 ex = imfilter(ex,filter);
%                 ey = imfilter(ey,filter);

                e_abs_av = e_abs_av + ex.^2 + ey.^2;
              
%                 
%                 e_abs_t(:,:,j) = e_abs(:, m*x_part_size+1:(m+1)*x_part_size);
                
            
            end
%             save([path4save, num2str(k,'%3.3d'), '_', num2str(m+1,'%1.1d'),'.dat'],'e_abs_t','-mat');

    
 res1 = e_abs_av/N;
 



