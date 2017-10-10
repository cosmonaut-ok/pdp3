function res1 = massive_fft;

dt = 1e-10;

w = 1.0e9;
Tmod = 2*pi/w;



path = 'G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\';
load('G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\ion_motion_16-48_10e12_01Ly_periodic.mat', 'param');
path4save = 'G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\_analisys\';


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
    
    t_start = 0.0e-7;
    t_end = 19.4e-7;
    dt = 1e-9;
    
    n_parts = 1;
    n_slices = 32;
    
    x_part_size = geometry.ngx/n_parts;

    e_abs_t = ones(geometry.ngy, geometry.ngx/n_parts, n_slices);
    rho_ions = ones(geometry.ngy, geometry.ngx/n_parts, n_slices);
    
    N = length(t_start:dt:t_end);
    
    
    k_max = floor(N/n_slices)

    
    for k = 10:k_max-1
        
        for m = 0:n_parts-1
            m
            j = 1;
            
            for t = (t_start + k*dt*n_slices):dt:(t_start + (k + 1)*dt*n_slices - dt)
                
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
                    
%                 rho = rho_el + rho_ions + rho_beam;    
%                 
%                 fi = field_3(rho,geometry,bc);
%                 [ex ey] = e_from_fi(fi,geometry,bc);
%                 e_abs = (ex.^2 + ey.^2);
%                 
%                 e_abs_t(:,:,j) = e_abs(:, m*x_part_size+1:(m+1)*x_part_size);
                 rho_ions(:,:,j) = rho_ions - rho_back;
                j = j + 1;
            
            end
%             save([path4save, num2str(k,'%3.3d'), '_', num2str(m+1,'%1.1d'),'.dat'],'e_abs_t','-mat');
             save([path4save, num2str(k,'%3.3d'), '_', num2str(m+1,'%1.1d'),'_ions','.dat'],'rho_ions','-mat');
            
            
        end
    end    
    
    
 res1 = 1;



