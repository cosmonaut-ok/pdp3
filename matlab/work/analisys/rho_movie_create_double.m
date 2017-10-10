function var = rho_movie_create_double;

w = 1.0e9;
Tmod = 2*pi/w;
h = ones(10)/100;

path_1 = 'G:\_results\mobile_16-48e14_05e11_01Ly_04.02\';
path_2 = 'G:\_results\moveless_16-48_10e12_01Ly_06_18\';

load('G:\_results\mobile_16-48e14_05e11_01Ly_04.02\mobile_16-48e14_05e11_01Ly.mat', 'param');
param_1 = param;

load('G:\_results\moveless_16-48_10e12_01Ly_06_18\moveless_16-48e14_10e12_01Ly.mat', 'param');
param_2 = param;

path4wr = 'G:\_results\moveless_16-48_10e12_01Ly_06_18\_movie\';

geometry_1 = param_1.geometry;
bc_1 = param_1.bc;
geometry_1 = calc_grid_step(geometry_1, bc_1);
rho_back_struct_1 = param_1.rho_back_struct;


geometry_2 = param_2.geometry;
bc_2 = param_2.bc;
geometry_2 = calc_grid_step(geometry_2, bc_2);
rho_back_struct_2 = param_2.rho_back_struct;



t_begin = 2.0e-7;
dt = 5e-10;
var = 'rho_sp';

% if rho_back_struct.enabled
    rho_back_1 = load_rho_back(rho_back_struct_1, geometry_1, bc_1)*1.6e-19;
    rho_back_2 = load_rho_back(rho_back_struct_2, geometry_2, bc_2)*1.6e-19;
% end

%-----determine visualization range-------

        clim1 = [0 10000];
        clim2 = [0 1000];
        clim3 = [-1e14 1e14];
        
        tick1 = [0 4000 8000];
        tick2 = [0 4000 8000];
        tick3 = [-5E13 0 5E13];
        
        tick_label_1 = ['   0'; '4000'; '8000'];
        tick_label_2 = ['  0'; '400'; '800'];
        tick_label_3 = ['-5';' 0'; ' 5'];
   
  
f = figure;
set(f, 'Position', [50 60 1050 700], 'Color', 'white');    
 
a1 = axes('Parent', f);
a2 = axes('Parent', f);
a3 = axes('Parent', f);
a4 = axes('Parent', f);

set(a1, 'Unit', 'normalized', 'Position', [0.06 0.73 0.8 0.2]);
set(a2, 'Unit', 'normalized', 'Position', [0.06 0.49 0.8 0.2]);
set(a3, 'Unit', 'normalized', 'Position', [0.06 0.25 0.8 0.2]);
set(a4, 'Unit', 'normalized', 'Position', [0.06 0.02 0.8 0.2]);

%initial---
z = zeros(geometry_1.ngy, geometry_1.ngx);
        
im1 = image(z,'Parent',a1, 'CDataMapping', 'scaled');
set(a1, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [' 0.0'; '0.15'; ' 0.3'; '0.45'; ' 0.6'],...
    'yTick', [128 256], 'yTickLabel', ['0.05'; ' 0.1'], 'Clim', clim1, 'YDir', 'normal'); 
im2 = image(z,'Parent',a2, 'CDataMapping', 'scaled');
set(a2, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim1, 'YDir', 'normal');
im3 = image(z,'Parent',a3, 'CDataMapping', 'scaled');
set(a3, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim2, 'YDir', 'normal');
im4 = image(z,'Parent',a4, 'CDataMapping', 'scaled');
set(a4, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim3, 'YDir', 'normal');
%---------------


%subscriptions----------
    set(get(a1,'XLabel'),'String','x, m', 'FontSize', 16.0);
    set(get(a1,'YLabel'),'String','y, m', 'FontSize', 16.0);

%     text_n0_0 = text(0,-0.5, 'n_i_,_e(x=0,t=0) = 3.2\cdot10^1^4 m^-^3', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_n0_L = text( 0,-0.8, 'n_i_,_e(x=L,t=0) = 3.2\cdot10^1^4 m^-^3', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_n0_L = text( 0,-1.1, 'T_e = 2.0 eV', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     
%     text_Te = text(0.3, -0.5, '\phi(x=0) = \phi(x = Lx) = 0; (Dirichlet bc)', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_Ti = text(0.3, -0.8, '\phi(y=0) = \phi(y = Ly) = 0 ; (Dirichlet bc)', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_Ti = text(0.3, -1.1, 'T_i = 0.0 eV', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     
%     text_f = text(0.7, -0.5, 'f_m_o_d = 1.6\cdot10^8 Hz', 'FontSize', 14, 'Parent', a1,'Unit','normalized'); 
%     text_n_b = text(0.7, -0.8, 'n_b_e_a_m = 1\cdot10^1^2 m^-^3', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_v_b = text(0.7, -1.1, 'v_b_e_a_m = 50\cdotv_T', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
% 
%     text_Nbe = text(0.95, -0.5, 'N_b_i_g_,_e_l = 1.2\cdot10^6', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_Nbi = text(0.95, -0.8, 'N_b_i_g_,_i_o_n = 0.0\cdot10^6', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
%     text_Ngr = text(0.95, -1.1, 'N_g_r_i_d = 2048\times256', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
    
    text_title = text(0.2, 1.2, strcat('Beam-plasma interaction,  t =','  ',num2str(0/Tmod,'%2.1f'),'  T_m_o_d'), ...
    'FontSize', 16, 'Parent', a1,'Unit','normalized');
    
%     text_n_el_1 = text(1.1, 0.9, '\Deltan_e_l,', 'Parent', a1,'Unit','normalized','FontSize', 16);
%     text_n_el_2 = text(1.1, 0.65, '10^1^4 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);
%     
%     text_n_ion_1 = text(1.1, 2.1, '\Deltan_i_o_n,', 'Parent', a1,'Unit','normalized','FontSize', 16);
%     text_n_ion_2 = text(1.1, 1.85, '10^1^4 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);
%     text_fi_1 = text(1.1, 2.1, '|E|_a_v_e_r,V/m', 'Parent', a1,'Unit','normalized','FontSize', 16);
%  
%     
%     text_n_beam_1 = text(1.1, 3.3, '\Deltan_b,', 'Parent', a1,'Unit','normalized','FontSize', 16);
%     text_n_beam_2 = text(1.1, 3.05, '10^1^2 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);

%   
    text1 = text(1.1, 0.7, '|E|, averaged', 'Parent', a1,'Unit','normalized');
    text2 = text(1.1, -0.5, '|E|, averaged', 'Parent', a1,'Unit','normalized');
    text3 = text(1.1, -1.7, '\nabla |E|', 'Parent', a1,'Unit','normalized');
    text4 = text(1.1, -2.9, '\Delta n_i', 'Parent', a1,'Unit','normalized');

    clrbr1 = colorbar('peer', a1,'Position', [0.88 0.73 0.03 0.2], 'Clim',clim1,'XTick',[], 'xTickLabel', [],...            %
       'YTick', tick1, 'yTickLabel', tick_label_1 );
    clrbr2 = colorbar('peer',a2, 'Position', [0.88 0.49 0.03 0.2],'Clim',clim2,'XTick',[], 'xTickLabel', [],...
        'YTick', tick1, 'yTickLabel', tick_label_1 );
    clrbr3 = colorbar('peer',a3, 'Position', [0.88 0.25 0.03 0.2],'Clim',clim3,'XTick',[], 'xTickLabel', [],...
        'YTick', tick2, 'yTickLabel', tick_label_2 );
    clrbr4 = colorbar('peer',a4, 'Position', [0.88 0.02 0.03 0.2],'Clim',clim3,'XTick',[], 'xTickLabel', [],...
        'YTick', tick3, 'yTickLabel', tick_label_3 );    
colormap('gray')

%----------


N = 13;
e_av_1 = zeros(geometry_1.ngy, geometry_1.ngx, N);
e_av_2 = zeros(geometry_2.ngy, geometry_2.ngx, N);

j = 1;

for k = 1:200

    tstart = (k-1)*15*dt + t_begin;
    tend = (k*15-1)*dt + t_begin;
    i = 1;

    for t = tstart:dt:tend
        
        %---------first loading-----------------------------
        
        load('-mat',strcat(path_1,'rho_electrons_',num2str(t),'.dat'),var) 
        rho_el = rho_sp;
        load('-mat',strcat(path_1,'rho_ions_',num2str(t),'.dat'),var) 
        rho_ions = rho_sp;     
        load('-mat',strcat(path_1,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;
        rho = rho_sp;
        rho = rho_el+rho_ions+rho_beam;
        
        fi = field_3(rho,geometry_1,bc_1);
        [ex ey] = e_from_fi(fi,geometry_1,bc_1);
        e_abs = (ex.^2 + ey.^2).^0.5;
        
        if mod(j,N) > 0
            e_av_1(:,:, mod(j,N)) = e_abs;
        else
            e_av_1(:,:, end) = e_abs;            
        end
        
        e_av = sum(e_av_1,3)/N;
        [grx gry] = gradient(e_av);
        
        abs_gr = (grx.^2 + gry.^2).^0.5;
        rho_ions = imfilter((rho_ions + rho_back_1*(-1))/1.6e-19, h);
        
        set(im2, 'CData', e_av);
        set(im3, 'CData', abs_gr);
        set(im4, 'CData', rho_ions);
        
        %----------------------------------------------------
        %---------second loading-----------------------------
        load('-mat',strcat(path_2,'rho_electrons_',num2str(t),'.dat'),var)    
        load('-mat',strcat(path_2,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;
        rho = rho_sp;
        rho = rho_el+rho_back_2+rho_beam;
        
        fi = field_3(rho,geometry_2,bc_2);
        [ex ey] = e_from_fi(fi,geometry_2,bc_2);
        e_abs = (ex.^2 + ey.^2).^0.5;
        
        if mod(j,N) > 0
            e_av_2(:,:, mod(j,N)) = e_abs;
        else
            e_av_2(:,:, end) = e_abs;            
        end
        
        e_av = sum(e_av_2,3)/N;
        [grx gry] = gradient(e_av);
        
        abs_gr = (grx.^2 + gry.^2).^0.5;
        
        set(im1, 'CData', e_av);
        
        j = j + 1;
        %-------------------------------------------------------

    set(text_title, 'String',strcat('Beam-plasma interaction,  t =',' ',num2str(t/Tmod,'%2.1f'),'  T_m_o_d', ';   1 - moveless, 2,3,4 - mobile ions'));
colormap('gray')

        figHandle = gcf;
        frame = getframe(figHandle);
        D(:,:,:,i)   = frame.cdata;
        D(:,:,:,i+1) = frame.cdata;
        D(:,:,:,i+2) = frame.cdata;
        i = i + 3;
   drawnow;
        
    end
    f2 = figure('Position', [20 20 100 100]);
    at = axes('Parent', f2);
    mov = immovie(D);
    movie2avi(mov, strcat(path4wr,'double_movie_',num2str(t/1e-7,'%3.2f'),'.avi'));
    clear D mov
    close(f2)
    pack



end
res = 1;
