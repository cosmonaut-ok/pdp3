function [res1 res2] = get_pulse_velocity;

% path = 'G:\_results\soliton_12e12\';
% load('G:\_results\soliton_12e12\soliton_22-42e14_Ly_10_12e12.mat', 'param');

% path = 'G:\_results\soliton_14e12\';
% load('G:\_results\soliton_14e12\soliton_22-42e14_Ly_10_14e12.mat', 'param');

path = 'G:\_results\_ion_motion_16-48_10e12_01Ly_periodic_05_12\';
load('G:\_results\_ion_motion_16-48_10e12_01Ly_periodic_05_12\ion_motion_16-48_10e12_01Ly_periodic.mat', 'param');

path = 'G:\_results\_ion_motion_16-48_14e12_01Ly_periodic_12_02\';
load('G:\_results\_ion_motion_16-48_14e12_01Ly_periodic_12_02\ion_motion_16-48_14e12_01Ly_periodic.mat', 'param');
% 
path = 'G:\_results\_ion_motion_16-48_18e12_01Ly_periodic_02_13\';
load('G:\_results\_ion_motion_16-48_18e12_01Ly_periodic_02_13\ion_motion_16-48_18e12_01Ly_periodic.mat', 'param');


filter = ones(9,12)/9/12;

geometry = param.geometry;
bc = param.bc;
geometry = calc_grid_step(geometry,bc);
rho_back_struct = param.rho_back_struct;

rho_back = load_rho_back(rho_back_struct, geometry, bc)*1.6e-19;

  j = 1;
        dt = 10e-9;
    tstart = 0.9e-6;
    tend = 1.4e-6;
    
    pos = 65;
%     pos = 670;
    
    res = length(tstart:dt:tend);
    
%     dt = 1e-9;
    
%         tstart = 0.27e-6*2;
%     tend = 0.27e-6*2 + 24*dt*0;

%     e_av = zeros(256,2048);
      rho_ions = zeros(256,1400);
%     e_x_av = zeros(256,2048);
%     e_y_av = zeros(256,2048);
    
    rho_ions = rho_back;
%     ne_av = zeros(256,2048);
%    e_max = ones(1, length(tstart:dt:tend));
    
    N = length(tstart:dt:tend)
%     fi_x_t = ones(N, geometry.ngx);
    rho_ions_y_t = ones(N, geometry.ngy);
%     rho_ions_x_t = ones(N, geometry.ngx);
%     e_abs_t = ones(geometry.ngy, geometry.ngx, N);
%     res1 = zeros(geometry.ngy,geometry.ngx);
%     e_abs_av = zeros(geometry.ngy,geometry.ngx);
%     res2 = res1;
%     res3 = res1;
    for t = tstart:dt:tend
        t 
        
        res1(j) = pos;
        try,
        load('-mat',strcat(path,'rho_electrons_',num2str(t),'.dat'),'rho_sp') 
        rho_el = rho_sp;


        load('-mat',strcat(path,'rho_ions_',num2str(t),'.dat'),'rho_sp') 
        rho_ions = rho_sp;     
        load('-mat',strcat(path,'rho_light_electrons_',num2str(t),'.dat'),'rho_sp') 
        rho_beam = rho_sp;
        catch,
            
            ttt = t+dt/20/10;
                load('-mat',strcat(path,'rho_electrons_',num2str(ttt),'.dat'),'rho_sp') 
        rho_el = rho_sp;

        load('-mat',strcat(path,'rho_ions_',num2str(ttt),'.dat'),'rho_sp') 
        rho_ions = rho_sp;     
        load('-mat',strcat(path,'rho_light_electrons_',num2str(ttt),'.dat'),'rho_sp') 
        rho_beam = rho_sp;
        end
%         rho = rho_el + rho_ions + rho_beam*0;

% 
% grid_position = round(27/60*2048);
% % %         
%         rho_ions_y_t(j,:) = (rho_ions(:,grid_position) - rho_back(:,grid_position))';

%         fi = field_5(rho,geometry,bc);
%         [ex ey] = e_from_fi(fi,geometry,bc);
%         
%         ex = imfilter(ex,filter);
        rho_ions = imfilter(rho_ions - rho_back,filter);


        
%         figure, plot(rho_ions(128,:))
%         figure, plot(rho_ions(:,945))
%         figure, imshow(rho_ions,[])
%         figure, imshow(rho_ions(:,900:1000),[])
%         afd


        
%         pos = near_left_maximum(rho_ions(128,:),pos);
        pos = near_left_maximum(rho_ions(:,945),pos);
%          pos = near_right_maximum(rho_ions(128,:),pos);
%         amplitude(j) = rho_ions(128,pos);
        amplitude(j) = rho_ions(pos,945);
        res1(j) = pos;
        
%         e_abs = (ex.^2 + ey.^2).^0.5;





        
%         e_x_av = e_x_av + ex;
%         e_y_av = e_y_av + ey;
%        e_abs_av = e_abs_av + e_abs;


   j = j+1;
    end

    res2 = amplitude;
    