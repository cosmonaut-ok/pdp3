function var = temp9(geometry,bc);

dt = 5e-10;

geometry = calc_grid_step(geometry,bc);

w = 1.0e9;
Tmod = 2*pi/w;

rho_back_struct(1) = struct('name', {'ions'},  ...
                            'charge', {1},     ...
                            'profile_fun',  {struct('handle',{@linear_fun}, 'param',{[2.4e14 4.0e14]})},...
                            'enabled', {1});

% h = ones(5,5)/25;

path = 'g:\_results\04.08\';
path4wr = 'g:\_results\04.08\_movie\';

t_begin = 6.75e-7;
var = 'rho_sp';


if rho_back_struct.enabled
    rho_back = load_rho_back(rho_back_struct, geometry, bc);
end



%-----determine visualization range-------


        clim1 = [-3.0e14 4.0e14];
        clim2 = [-4.0e14 4.0e14];
        clim3 = [-0.5e4 0];
   
  
f = figure;
set(f, 'Position', [50 60 1050 700], 'Color', 'white');    
 
a1 = axes('Parent', f);
a2 = axes('Parent', f);
a3 = axes('Parent', f);
  set(a1, 'Unit', 'normalized', 'Position', [0.06 0.25 0.8 0.2]);
    set(a2, 'Unit', 'normalized', 'Position', [0.06 0.49 0.8 0.2]);
    set(a3, 'Unit', 'normalized', 'Position', [0.06 0.73 0.8 0.2]);
    
%initial---
t = 0;
        load('-mat',strcat(path,'rho_electrons_',num2str(t),'.dat'),var) 
        rho_el = rho_sp;
        load('-mat',strcat(path,'rho_ions_',num2str(t),'.dat'),var) 
        rho_ions = rho_sp;     
        load('-mat',strcat(path,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;

        
        rho_el = rho_el/1.6e-19 + rho_back;
        rho_ions = rho_ions/1.6e-19 + rho_back*(-1);
        rho_beam = rho_beam/1.6e-19;
        
        im1 = image(rho_el,'Parent',a1, 'CDataMapping', 'scaled');
        set(a1, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [' 0.0'; '0.15'; ' 0.3'; '0.45'; ' 0.6'],...
                'yTick', [128 256], 'yTickLabel', ['0.05'; ' 0.1'], 'Clim', clim1, 'YDir', 'normal'); 
        im2 = image(rho_ions,'Parent',a2, 'CDataMapping', 'scaled');
        set(a2, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim2, 'YDir', 'normal');
        im3 = image(rho_beam,'Parent',a3, 'CDataMapping', 'scaled');
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
    text_n_el_1 = text(1.1, 0.9, '\Deltan_e_l,', 'Parent', a1,'Unit','normalized','FontSize', 16);
    text_n_el_2 = text(1.1, 0.65, '10^1^4 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);
    
    text_n_ion_1 = text(1.1, 2.1, '\Deltan_i_o_n,', 'Parent', a1,'Unit','normalized','FontSize', 16);
    text_n_ion_2 = text(1.1, 1.85, '10^1^4 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);
    
    text_n_beam_1 = text(1.1, 3.3, '|E|', 'Parent', a1,'Unit','normalized','FontSize', 16);
    text_n_beam_2 = text(1.1, 3.05, '10^3 V/m', 'Parent', a1,'Unit','normalized','FontSize', 12);

%     
    clrbr1 = colorbar('peer', a1,'Position', [0.88 0.25 0.03 0.2], 'Clim',clim1,'XTick',[], 'xTickLabel', [],...            %
       'YTick',[-2.0E14 0 2.0E14], 'yTickLabel', ['-2.0'; ' 0.0'; ' 2.0'] );
    clrbr2 = colorbar('peer',a2, 'Position', [0.88 0.49 0.03 0.2],'Clim',clim2,'XTick',[], 'xTickLabel', [],...
        'YTick',[-2.0E14 0 2.0E14], 'yTickLabel', ['-2.0'; ' 0.0'; ' 2.0'] );
    clrbr3 = colorbar('peer',a3, 'Position', [0.88 0.73 0.03 0.2],'Clim',clim3,'XTick',[], 'xTickLabel', [],...
        'YTick',[-1E3 -3E3], 'yTickLabel', ['3'; '1'] );
    
colormap('gray')

%----------
h = ones(10)/100;

for k = 1:100

    tstart = (k-1)*25*dt + t_begin;
    tend = (k*25-1)*dt + t_begin;
    i = 1;

    for t = tstart:dt:tend
        load('-mat',strcat(path,'rho_electrons_',num2str(t),'.dat'),var) 
        rho_el = rho_sp;
        load('-mat',strcat(path,'rho_ions_',num2str(t),'.dat'),var) 
        rho_ions = rho_sp;     
        load('-mat',strcat(path,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;
        rho = rho_sp;
        rho = rho_el+rho_ions+rho_beam;
        
        rho_el = rho_el/1.6e-19 + rho_back;
        rho_el = rho_el*(-1);
        rho_ions = rho_ions/1.6e-19 + rho_back*(-1);
        rho_beam = rho_beam/1.6e-19;
        
        fi = field_3(rho,geometry,bc);
        [ex ey] = e_from_fi(fi,geometry,bc);
        e_abs = (ex.^2 + ey.^2).^0.5;
        
        rho_el = imfilter(rho_el,h);
        rho_ions = imfilter(rho_ions,h);
        
        set(im1, 'CData', rho_el);
        set(im2, 'CData', rho_ions);
        set(im3, 'CData', (-1)*e_abs);
    set(text_title, 'String',strcat('Beam-plasma interaction,  t =','  ',num2str(t/Tmod,'%2.1f'),'  T_m_o_d'));
adfafd
    
colormap('gray')

        figHandle = gcf;
        frame = getframe(figHandle);
        D(:,:,:,i)   = frame.cdata;
        D(:,:,:,i+1) = frame.cdata;
        D(:,:,:,i+2) = frame.cdata;
        i = i + 3;
        
    end
    f2 = figure('Position', [20 20 100 100]);
    at = axes('Parent', f2);
    mov = immovie(D);
    movie2avi(mov, strcat(path4wr,'rho_movie_',num2str(t/1e-7,'%3.2f'),'.avi'));
    clear D mov
    close(f2)
    pack



end
res = 1;
