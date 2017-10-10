function var = rho_movie_create;

w = 1.0e9;
Tmod = 2*pi/w;

filter = ones(9,12)/9/12;

% path = 'D:\_results\mobile_1.6e14_10
% path = 'G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\';
% path4wr = 'G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\_movie\';
% load('G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\ion_motion_16-48_10e12_01Ly_periodic.mat', 'param');
%
% path = 'G:\_results\ion_motion_20-44e14_10e12_01Ly_periodic_20_12\';
% path4wr = 'G:\_results\ion_motion_20-44e14_10e12_01Ly_periodic_20_12\_movie\';
% load('G:\_results\ion_motion_20-44e14_10e12_01Ly_periodic_20_12\ion_motion_20-44e14_10e12_01Ly_periodic.mat', 'param');

% path = 'G:\_results\ion_motion_16-48e14_2e11_01Ly_periodic_29_12\';
% path4wr = 'G:\_results\ion_motion_16-48e14_2e11_01Ly_periodic_29_12\_movie\';
% load('G:\_results\ion_motion_16-48e14_2e11_01Ly_periodic_29_12\ion_motion_16-48e14_2e11_01Ly_periodic.mat', 'param');

% path = 'I:\_results\ion_motion_16-48_2e11_10Ly_periodic_06_03\';
% path4wr = 'I:\_results\ion_motion_16-48_2e11_10Ly_periodic_06_03\_movie\';
% load('I:\_results\ion_motion_16-48_2e11_10Ly_periodic_06_03\ion_motion_16-48_2e11_10Ly_periodic.mat', 'param');
% 
% path = 'I:\_results\ion_motion_16-48_6e11_10Ly_periodic_11_03\';
% path4wr = 'I:\_results\ion_motion_16-48_6e11_10Ly_periodic_11_03\\_movie\';
% load('I:\_results\ion_motion_16-48_6e11_10Ly_periodic_11_03\ion_motion_16-48_6e11_10Ly_periodic.mat', 'param');
% 
% path = 'I:\_results\ion_motion_16-48_6e10_10Ly_periodic_14_03\';
% path4wr = 'I:\_results\ion_motion_16-48_6e10_10Ly_periodic_14_03\_movie\';
% load('I:\_results\ion_motion_16-48_6e10_10Ly_periodic_14_03\ion_motion_16-48_6e10_10Ly_periodic.mat', 'param');

path = 'D:\masha\pdp_2\';
path4wr = 'D:\masha\pdp_2\movie\';
load('D:\masha\pdp_2\config.mat', 'param');

geometry = param.geometry;
bc = param.bc;

geometry = calc_grid_step(geometry, bc);
rho_back_struct = param.rho_back_struct;

t_begin = 1e-10;
dt = 1e-10;
var = 'rho_beam_electrons';

% if rho_back_struct.enabled
rho_back = load_rho_back(rho_back_struct, geometry, bc)*1.6e-19;

% end0

%-----determine visualization range-------



clim1 = [-3e14 3e14];
clim2 = [-3e14 3e14];
clim3 = [-0.3e13 0.3e13];
clim3 = [-0.6e12 0.6e12];

tick1 = [-2.0E14 0 2.0E14];
tick2 = [-2.0E14 0 2.0E14];
tick3 = [-2E12 0 2E12];

tick_label_1 = ['-2'; ' 0'; ' 2'];
tick_label_2 = ['-2'; ' 0'; ' 2'];
tick_label_3 = ['-2';' 0'; ' 2'];



clim1 = [-1e14 0];
tick1 = [0 1e14 2e14 3e14];
tick_label_1 = ['0'; '1'; '2' ; '3'];


f = figure;
set(f, 'Position', [50 60 1050 700], 'Color', 'white');

a1 = axes('Parent', f);
a2 = axes('Parent', f);
a3 = axes('Parent', f);
% a4 = axes('Parent', f);

set(a1, 'Unit', 'normalized', 'Position', [0.06 0.25 0.8 0.2]);
set(a2, 'Unit', 'normalized', 'Position', [0.06 0.49 0.8 0.2]);
set(a3, 'Unit', 'normalized', 'Position', [0.06 0.73 0.8 0.2]);
% set(a4, 'Unit', 'normalized', 'Position', [0.06 0.02 0.8 0.2]);

%initial---
z = zeros(geometry.ngy, geometry.ngx);

im1 = image(z,'Parent',a1, 'CDataMapping', 'scaled');
set(a1, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [' 0.0'; '0.15'; ' 0.3'; '0.45'; ' 0.6'],...
    'yTick', [128 256], 'yTickLabel', ['0.05'; ' 0.1'], 'Clim', clim1, 'YDir', 'normal');
im2 = image(z,'Parent',a2, 'CDataMapping', 'scaled');
set(a2, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim2, 'YDir', 'normal');
im3 = image(z,'Parent',a3, 'CDataMapping', 'scaled');
set(a3, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim3, 'YDir', 'normal');
%---------------


%subscriptions----------
set(get(a1,'XLabel'),'String','x, m', 'FontSize', 16.0);
set(get(a1,'YLabel'),'String','y, m', 'FontSize', 16.0);

text_n0_0 = text(0,-0.5, 'n_i_,_e(x=0,t=0) = 1.6\cdot10^1^4 m^-^3', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
text_n0_L = text( 0,-0.8, 'n_i_,_e(x=L,t=0) = 4.8\cdot10^1^4 m^-^3', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
text_n0_L = text( 0,-1.1, 'T_e = 2.0 eV', 'FontSize', 14, 'Parent', a1,'Unit','normalized');

text_Te = text(0.3, -0.5, '\phi(x=0) = \phi(x = Lx) = 0; (Dirichlet bc)', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
text_Ti = text(0.3, -0.8, '\phi(y=0) = \phi(y = Ly) ; (periodic bc)', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
text_Ti = text(0.3, -1.1, 'T_i = 0.2 eV', 'FontSize', 14, 'Parent', a1,'Unit','normalized');

text_f = text(0.7, -0.5, 'f_m_o_d = 1.6\cdot10^8 Hz', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
text_n_b = text(0.7, -0.8, 'n_b_e_a_m = 6\cdot10^1^0 m^-^3', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
text_v_b = text(0.7, -1.1, 'v_b_e_a_m = 50\cdotv_T', 'FontSize', 14, 'Parent', a1,'Unit','normalized');

text_Nbe = text(0.95, -0.5, 'N_b_i_g_,_e_l = 2.0\cdot10^6', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
text_Nbi = text(0.95, -0.8, 'N_b_i_g_,_i_o_n = 2.0\cdot10^6', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
text_Ngr = text(0.95, -1.1, 'N_g_r_i_d = 2048\times256', 'FontSize', 14, 'Parent', a1,'Unit','normalized');

text_title = text(0.4, 3.65, strcat('Beam-plasma interaction,  \omega_m_o_d\cdott =','  ',num2str(0/Tmod,'%2.1f')), ...
    'FontSize', 16, 'Parent', a1,'Unit','normalized');



text_n_ion_1 = text(1.07, 2.1, '\Deltan_i_o_n', 'Parent', a1,'Unit','normalized','FontSize', 24);
%     text_n_ion_1 = text(1.07, 2.1, 'n_i_o_n', 'Parent', a1,'Unit','normalized','FontSize', 24);
text_n_ion_2 = text(1.07, 1.65, '10^1^4 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);

%     text_n_el_1 = text(1.07, 0.9, '\Deltan_e_l', 'Parent', a1,'Unit','normalized','FontSize', 24);
%     text_n_el_2 = text(1.07, 0.45, '10^1^4 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);

text_n_el_1 = text(1.1, 0.9, 'E_a_b_s', 'Parent', a1,'Unit','normalized','FontSize', 24);
text_n_el_2 = text(1.1, 0.45, '10^3 V/m', 'Parent', a1,'Unit','normalized','FontSize', 12);
%
%     text_n_ion_1 = text(1.1, 2.1, '\phi', 'Parent', a1,'Unit','normalized','FontSize', 24);
%     text_n_ion_2 = text(1.1, 1.65, 'V', 'Parent', a1,'Unit','normalized','FontSize', 12);



text_n_beam_1 = text(1.07, 3.3, 'n_b_e_a_m', 'Parent', a1,'Unit','normalized','FontSize', 24);
text_n_beam_2 = text(1.07, 2.85, '10^1^2 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);

%
clrbr1 = colorbar('peer', a1,'Position', [0.88 0.25 0.01 0.2], 'Clim',clim1,'XTick',[], 'xTickLabel', [],...            %
    'YTick', tick1, 'yTickLabel', tick_label_1 );
clrbr2 = colorbar('peer',a2, 'Position', [0.88 0.49 0.01 0.2],'Clim',clim2,'XTick',[], 'xTickLabel', [],...
    'YTick', tick2, 'yTickLabel', tick_label_2 );
clrbr3 = colorbar('peer',a3, 'Position', [0.88 0.73 0.01 0.2],'Clim',clim3,'XTick',[], 'xTickLabel', [],...
    'YTick', tick3, 'yTickLabel', tick_label_3 );

colormap('gray')

%----------


% N = 25;
% ex_av = zeros(geometry.ngy, geometry.ngx, N);
% ey_av = zeros(geometry.ngy, geometry.ngx, N);
j = 1;

N = 20;




for k = 1:10

    tstart = (k-1)*N*dt + t_begin;
    tend = (k*N-1)*dt + t_begin;
    i = 1;

    for t = tstart:dt:tend

        try,
            load('-mat',strcat(path,'rho_beam_electrons_',num2str(t),'.dat'),'rho_sp')
            
            %         load('-mat',strcat(path,'rho_electrons_',num2str(t),'.dat'),'rho_sp')
            %         rho_el = rho_sp;
            %
            %         load('-mat',strcat(path,'rho_ions_',num2str(t),'.dat'),var)
            %         rho_ions = rho_sp;
            %         load('-mat',strcat(path,'rho_light_electrons_',num2str(t),'.dat'),var)
            %         rho_beam = rho_sp;
        catch,

            ttt = t+dt/20;
            load('-mat',strcat(path,'rho_beam_electrons_',num2str(ttt),'.dat'),'rho_sp')
            %             load('-mat',strcat(path,'rho_electrons_',num2str(ttt),'.dat'),'rho_sp')
            %             rho_el = rho_sp;
            %
            %             load('-mat',strcat(path,'rho_ions_',num2str(ttt),'.dat'),var)
            %             rho_ions = rho_sp;
            %             load('-mat',strcat(path,'rho_light_electrons_',num2str(ttt),'.dat'),var)
            %             rho_beam = rho_sp;
        end
      min(min(rho_sp))
      max(max(rho_sp))

        %         rho = rho_sp;
%         rho = rho_el+rho_ions+rho_beam;
        %         rho = rho_el + rho_beam + rho_back;


        %         rho_el = (rho_el + rho_back)/1.6e-19;
        %         rho_el = rho_el*(-1);
        %         rho_ions = rho_ions/1.6e-19 + rho_back*(-1);
%         rho_beam = rho_beam/1.6e-19;

%         fi = field_3(rho,geometry,bc);
%         [ex ey] = e_from_fi(fi,geometry,bc);
%         e_abs = (ex.^2 + ey.^2).^0.5;



        %         if mod(j,N) > 0
        %             ex_av(:,:, mod(j,N)) = imfilter(ex,filter);
        %             ey_av(:,:, mod(j,N)) = imfilter(ey,filter);
        % %             rho_av(:,:, mod(j,N)) = ((-1)*rho_el - rho_back)/1.6e-19;
        %         else
        %             ex_av(:,:, end) = imfilter(ex,filter);
        %             ey_av(:,:, end) = imfilter(ey,filter);
        % %             e_av(:,:, end) = e_abs;
        % %             rho_av(:,:, end) = ((-1)*rho_el - rho_back)/1.6e-19;
        %         end
        j = j + 1;



        %         set(im1, 'CData', sum(ex_av,3)/N);
        %         set(im2, 'CData', sum(ey_av,3)/N);
        %         set(im3, 'CData', ((sum(ex_av,3)).^2 + (sum(ex_av,3)).^2).^0.5/N);

%         rho_el_delta = (rho_el + rho_back)/1.6e-19*(-1);
%         rho_ion_delta = (rho_ions - rho_back)/1.6e-19;

        %
        %         set(im1, 'CData', imfilter(rho_el_delta,filter));
        set(im1, 'CData', rho_sp/1.6e-19);

%         set(im2, 'CData', imfilter(rho_ion_delta,filter));

        %         set(im2, 'CData', fi);
%         set(im3, 'CData', imfilter(rho_beam,filter)*(-1));


        set(text_title, 'String',strcat('Beam-plasma interaction, non-isothermal plasma  \omega_0\cdott =  ','  ',num2str(t/Tmod,'%2.1f')));
        %         p = plot(sum(rho_beam)/256*10*(-1), 'Parent', a1);
        %         set(a1, 'ylim', [0 2e13], 'xlim', [0 2048])
        %         p = plot(fi(128,:), 'Parent', a4);
        %         set(a4, 'ylim', [-150 150], 'xlim', [0 2048])

%         colormap('jet')



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
    movie2avi(mov, strcat(path4wr,'field_movie_',num2str(t/1e-7,'%3.2f'),'.avi'));
    clear mov
    close(f2)
    pack



end
var = 1;
