function [res1 res2] = show_single_fig;

dt = 1e-10;

w = 1.0e9;
Tmod = 2*pi/w;

filter = ones(9,12)/9/12;
% filter = ones(3,4)/12;

% path = 'G:\_results\mobile_48-16e14_10e12_01Ly_10_12\';
% load('G:\_results\mobile_48-16e14_10e12_01Ly_10_12\mobile_48-16e14_10e12_01Ly.mat', 'param');

% path = 'G:\_results\bd_mobile_32e14_10e12_01Ly_periodic_critical_09_14\';
% load('G:\_results\bd_mobile_32e14_10e12_01Ly_periodic_critical_09_14\bd_mobile_32e14_10e12_01Ly_periodic_critical.mat', 'param');

% path = 'G:\_results\bd_48e14_10e12_01Ly_periodic_10_09_supercrit\';
% load('G:\_results\bd_48e14_10e12_01Ly_periodic_10_09_supercrit\bd_48e14_10e12_01Ly_periodic_supercrit.mat', 'param');

% path = 'G:\_results\bd_mobile_16e14_10e12_01Ly_periodic_subcrit_08_31\';
% load('G:\_results\bd_mobile_16e14_10e12_01Ly_periodic_subcrit_08_31\bd_mobile_16e14_10e12_01Ly_periodic_subcrit.mat', 'param');

% path = 'G:\_results\mobile_20-44e14_10e12_01Ly_08_09\';
% load('G:\_results\mobile_20-44e14_10e12_01Ly_08_09\mobile_20-44e14_10e12_01Ly_precise.mat', 'param');

% path = 'G:\_results\_mobile_16-48e14_10e12_01Ly_07_19\';
% load('G:\_results\_mobile_16-48e14_10e12_01Ly_07_19\mobile_16-48e14_10e12_01Ly.mat', 'param');

% path = 'G:\_results\bd_48e14_10e12_01Ly_periodic_30_11_supercrit\bd_48e14_10e12_01Ly_periodic_30_11_supercrit\';
% load('G:\_results\bd_48e14_10e12_01Ly_periodic_30_11_supercrit\bd_48e14_10e12_01Ly_periodic_30_11_supercrit\bd_48e14_10e12_01Ly_periodic_supercrit.mat', 'param');

% 
% path = 'G:\_results\moveless_16-48_03e11_10Ly_06_20\';
% load('G:\_results\moveless_16-48_03e11_10Ly_06_20\moveless_16-48e14_03e11_01Ly.mat', 'param');

% path = 'G:\_results\moveless_16-48e14_10e12_1Ly_06_16\';
% load('G:\_results\moveless_16-48e14_10e12_1Ly_06_16\moveless_16-48e14_10e12_1Ly.mat', 'param');

% path = 'G:\_results\moveless_16-48e14_10e12_01Ly_06_22\';
% load('G:\_results\moveless_16-48e14_10e12_01Ly_06_22\moveless_16-48e14_10e12_01Ly.mat', 'param');

% path = 'G:\_results\moveless_16-48_10e12_01Ly_06_18\';
% load('G:\_results\moveless_16-48_10e12_01Ly_06_18\moveless_16-48e14_10e12_01Ly.mat', 'param');

% path = 'G:\_results\mobile_16-48e14_10e12_1Ly_06_24\';
% load('G:\_results\mobile_16-48e14_10e12_1Ly_06_24\mobile_16-48e14_10e12_1Ly.mat', 'param');

% path = 'G:\_results\mobile_16-48e14_10e12_01Ly_07_09_precise\';
% load('G:\_results\mobile_16-48e14_10e12_01Ly_07_09_precise\mobile_16-48e14_10e12_01Ly.mat', 'param');

% path = 'G:\_results\moveless_16-48e14_01e11_10Ly_07_23\';
% load('G:\_results\moveless_16-48e14_01e11_10Ly_07_23\moveless_16-48e14_01e11_10Ly', 'param');

% path = 'G:\_results\moveless_16-48e14_05e10_10Ly_07_24\';
% load('G:\_results\moveless_16-48e14_05e10_10Ly_07_24\moveless_16-48e14_05e10_10Ly.mat', 'param');

% path = 'G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\';
% load('G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\ion_motion_16-48_10e12_01Ly_periodic.mat', 'param');

% path = 'G:\_results\mobile_16-48e14_10e12_01_Ly_12_01_08\';
% load('G:\_results\mobile_16-48e14_10e12_01_Ly_12_01_08\mobile_16-48e14_10e12_01_Ly_03_Ti.mat', 'param');

% path = 'I:\_results\moveless_16-48e14_wide_beam10e12_07_17\';
% load('I:\_results\moveless_16-48e14_wide_beam10e12_07_17\moveless_16-48e14_wide_beam10e12_07_17.mat', 'param');
% 
% path = 'I:\_results\moveless_16-48e14_thin_beam10e12_07_21\';
% load('I:\_results\moveless_16-48e14_thin_beam10e12_07_21\moveless_16-48e14_thin_beam10e12_07_21.mat', 'param');
% 
% 
% 
% path = 'I:\_results\moveless_48-16e14_10e12_01Ly_07_27\';
% load('I:\_results\moveless_48-16e14_10e12_01Ly_07_27\moveless_48-16e14_thin_beam10e12_07_27.mat', 'param');
% 
% path = 'I:\_results\disser_mobile_16-48_1e12_thin_07_29\';
% load('I:\_results\disser_mobile_16-48_1e12_thin_07_29\ion_motion_16-48_10e12_01Ly_periodic.mat', 'param');


% path = 'I:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\';
% load('I:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\ion_motion_16-48_10e12_01Ly_periodic.mat', 'param');


% path = 'G:\_results\_ion_motion_16-48_10e12_01Ly_periodic_05_12\';
% load('G:\_results\_ion_motion_16-48_10e12_01Ly_periodic_05_12\ion_motion_16-48_10e12_01Ly_periodic.mat', 'param');

% path = 'G:\_results\_ion_motion_16-48_14e12_01Ly_periodic_12_02\';
% load('G:\_results\_ion_motion_16-48_14e12_01Ly_periodic_12_02\ion_motion_16-48_14e12_01Ly_periodic.mat', 'param');

path = 'G:\_results\_ion_motion_16-48_18e12_01Ly_periodic_02_13\';
load('G:\_results\_ion_motion_16-48_18e12_01Ly_periodic_02_13\ion_motion_16-48_18e12_01Ly_periodic.mat', 'param');


geometry = param.geometry;
bc = param.bc;
geometry = calc_grid_step(geometry,bc);
rho_back_struct = param.rho_back_struct;


t_begin = 2.0e-8;
var = 'rho_sp';
var1 = 'rho';


% if rho_back_struct.enabled
    rho_back = load_rho_back(rho_back_struct, geometry, bc)*1.6e-19;
% end



%-----determine visualization range-------


        clim2 = [-1.0e13 0];
        clim3 = [-4.0e14 4.0e14];
        clim1 = [-1e4 1e4];
   
  
f = figure;
set(f, 'Position', [50 60 1050 525], 'Color', 'white');    
 
a1 = axes('Parent', f);
a2 = axes('Parent', f);
a3 = axes('Parent', f);
  set(a1, 'Unit', 'normalized', 'Position', [0.06 0.15 0.8 0.25]);
    set(a2, 'Unit', 'normalized', 'Position', [0.06 0.43 0.8 0.25]);
    set(a3, 'Unit', 'normalized', 'Position', [0.06 0.72 0.8 0.25]);
    
%initial---

        z = zeros(geometry.ngy, geometry.ngx);
        
        im1 = image(z,'Parent',a1, 'CDataMapping', 'scaled');
        set(a1, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [' 0.0'; '0.15'; ' 0.3'; '0.45'; ' 0.6'],...
                'yTick', [128 256], 'yTickLabel', ['0.05'; ' 0.1'], 'Clim', clim1, 'YDir', 'normal', 'FontSize', 16); 
        im2 = image(z,'Parent',a2, 'CDataMapping', 'scaled');
        set(a2, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim2, 'YDir', 'normal');
        im3 = image(z,'Parent',a3, 'CDataMapping', 'scaled');
        set(a3, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim3, 'YDir', 'normal');



        
        %----------------------------------
%     set(get(a1,'XLabel'),'String','x, m', 'FontSize', 16.0);
%     set(get(a1,'YLabel'),'String','y, m', 'FontSize', 16.0);
% 
%     text_n0_0 = text(0,-0.5, 'n_i_,_e(x=0,t=0) = 2.4\cdot10^1^4 m^-^3', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_n0_L = text( 0,-0.8, 'n_i_,_e(x=L,t=0) = 4.0\cdot10^1^4 m^-^3', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_n0_L = text( 0,-1.1, 'T_e = 2.0 eV', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     
%     text_Te = text(0.3, -0.5, '\phi(x=0) = \phi(x=Lx) = 0; (Dirichlet bc)', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_Ti = text(0.3, -0.8, '\phi(y=0) = \phi(y = Ly); (periodic bc)', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_Ti = text(0.3, -1.1, 'T_i = 0.2 eV', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     
%     text_f = text(0.7, -0.5, 'f_m_o_d = 1.6\cdot10^8 Hz', 'FontSize', 14, 'Parent', a1,'Unit','normalized'); 
%     text_n_b = text(0.7, -0.8, 'n_b_e_a_m = 5\cdot10^1^1 m^-^3', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_v_b = text(0.7, -1.1, 'v_b_e_a_m = 50\cdotv_T', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
% 
%     text_Nbe = text(0.95, -0.5, 'N_b_i_g_,_e_l = 1.2\cdot10^6', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_Nbi = text(0.95, -0.8, 'N_b_i_g_,_i_o_n = 4.0\cdot10^5', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_Ngr = text(0.95, -1.1, 'N_g_r_i_d = 2048\times256', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     
%     text_title = text(0.4, 3.65, strcat('Beam-plasma interaction,  t =','  ',num2str(t/Tmod,'%2.1f'),'  T_m_o_d'), ...
%     'FontSize', 16, 'Parent', a1,'Unit','normalized');
%     
    text_n_el_1 = text(1.1, 0.9, '\phi, V', 'Parent', a1,'Unit','normalized','FontSize', 16);
%     text_n_el_2 = text(1.1, 0.65, '10^1^4 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);
    
    text_n_ion_1 = text(1.1, 2.1, 'n_b_e_a_m,', 'Parent', a1,'Unit','normalized','FontSize', 16);
    text_n_ion_2 = text(1.1, 1.85, '10^1^2 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);
    
    text_n_beam_1 = text(1.1, 3.3, '|E|', 'Parent', a1,'Unit','normalized','FontSize', 16);
    text_n_beam_2 = text(1.1, 3.05, '10^3 V/m', 'Parent', a1,'Unit','normalized','FontSize', 12);

%     
    clrbr1 = colorbar('peer', a1,'Position', [0.88 0.15 0.03 0.25], 'Clim',clim1,'XTick',[], 'xTickLabel', [],...            %
       'YTick',[-20 0 20], 'yTickLabel', ['-20'; '  0'; ' 20'], 'FontSize', 16 );
    clrbr2 = colorbar('peer',a2, 'Position', [0.88 0.43 0.03 0.25],'Clim',clim2,'XTick',[], 'xTickLabel', [],...
        'YTick',[-8E12 -4E12], 'yTickLabel', ['-8'; '-4'], 'FontSize', 16 );
    clrbr3 = colorbar('peer',a3, 'Position', [0.88 0.72 0.03 0.25],'Clim',clim3,'XTick',[], 'xTickLabel', [],...
        'YTick',[-1E3 -3E3], 'yTickLabel', ['3'; '1'] );
    
colormap('gray')

%----------
    
    j = 1;
        dt = 1e-9;
    tstart = 1.193e-6;
    tend = 1.193e-6 + 24*dt*0;
%     dt = 1e-9;
    
%         tstart = 0.27e-6*2;
%     tend = 0.27e-6*2 + 24*dt*0;

%     e_av = zeros(256,2048);
    e_x_av = zeros(256,2048);
    e_y_av = zeros(256,2048);
    
    rho_ions = rho_back;
%     ne_av = zeros(256,2048);
%    e_max = ones(1, length(tstart:dt:tend));
    
    N = length(tstart:dt:tend)
%     fi_x_t = ones(N, geometry.ngx);
    rho_ions_y_t = ones(N, geometry.ngy);
%     rho_ions_x_t = ones(N, geometry.ngx);
%     e_abs_t = ones(geometry.ngy, geometry.ngx, N);
    res1 = zeros(geometry.ngy,geometry.ngx);
    e_abs_av = zeros(geometry.ngy,geometry.ngx);
    res2 = res1;
%     res3 = res1;
    for t = tstart:dt:tend
        t 
%         try,
        load('-mat',strcat(path,'rho_electrons_',num2str(t),'.dat'),'rho_sp') 
        rho_el = rho_sp;

        load('-mat',strcat(path,'rho_ions_',num2str(t),'.dat'),var) 
        rho_ions = rho_sp;     
        load('-mat',strcat(path,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;
%         catch,
            
%             ttt = t+dt/20;
%                 load('-mat',strcat(path,'rho_electrons_',num2str(ttt),'.dat'),'rho_sp') 
%         rho_el = rho_sp;
% 
%         load('-mat',strcat(path,'rho_ions_',num2str(ttt),'.dat'),var) 
%         rho_ions = rho_sp;     
%         load('-mat',strcat(path,'rho_light_electrons_',num2str(ttt),'.dat'),var) 
%         rho_beam = rho_sp;
%         end
        rho = rho_el + rho_ions + rho_beam;

% 
% grid_position = round(27/60*2048);
% % %         
%         rho_ions_y_t(j,:) = (rho_ions(:,grid_position) - rho_back(:,grid_position))';

%         fi = field_5(rho,geometry,bc);
%         [ex ey] = e_from_fi(fi,geometry,bc);
%         e_abs = (ex.^2 + ey.^2).^0.5;
%          e_abs = (ex.^2 + ey.^2).^0.5;
%          e_abs = imfilter(e_abs,filter);
        
%         e_x_av = e_x_av + ex;
%         e_y_av = e_y_av + ey;
res1 = res1 + imfilter(rho_ions,filter) - rho_back;
%        e_abs_av = e_abs_av + e_abs;


   j = j+1;
    end
    
%     res1 = imfilter(rho_ions,filter) - rho_back*1;
% res1 = e_abs_av/N;
%     res2 = 1;

% res1 = e_x_av/N;
% res2 = e_y_av/N;

%     res1 = imfilter(rho_ions,filter) - rho_back*1;
% res1 = e_abs_av/N;
% res1 = e_x_av;
% res2 = e_y_av;
% 
% res1 = rho_ions - rho_back;

res2 = res1;

end

