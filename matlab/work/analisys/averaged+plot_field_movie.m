function res = averaged+plot_field_movie
dt = 5e-10;
%


geometry.x_size = 0.6;
geometry.y_size = 0.1;
geometry.ngx = 2048;
geometry.ngy = 256;

geometry.dx = geometry.x_size/(geometry.ngx - 1);
geometry.dy = geometry.y_size/geometry.ngy;

w = 1.0e9;
Tmod = 2*pi/w;

bc.x_type = 'dirichlet';
bc.y_type = 'periodic';

rho_back_struct(1) = struct('name', {'ions'},  ...
                            'charge', {1},     ...
                            'profile_fun',  {struct('handle',{@linear_fun}, 'param',{[2.4e14 4.0e14]})},...
                            'enabled', {1});

% xSize = 2e-0;
% ySize = 0.5e-0;
% 
% ngx = 512;
% ngy = 128;

h = ones(5,5)/25;

path = 'd:\_results\03.06\';
path4wr = 'd:\_results\03.06\movie\';
% path = 'c:\movie';
var = 'rho_sp';
% name_sp = 'rho_beam_';

t_begin = 0e-9;

% 
% tit = title(strcat('t =  ',num2str(t_begin/1e-9,'%4.1f'), ', ns'));
% set(tit, 'FontSize', 14);

if rho_back_struct.enabled
    rho_back = load_rho_back(rho_back_struct, geometry, bc)*1.6e-19;
end




        clim1 = [0 13e3];
        clim2 = [0 13e3];
        clim3 = [0 120];
        
        Tick_position_1 = [4000 8000 12000];
        Tick_position_2 = [4000 8000 12000];
        
        Tick_label_1 = ['1000'; '3000'; '5000'];
        Tick_label_2 = ['1000'; '3000'; '5000'];
   
  
f = figure;
set(f, 'Position', [50 60 1050 700], 'Color', 'white');    
 
a1 = axes('Parent', f);
a2 = axes('Parent', f);
a3 = axes('Parent', f);
  set(a1, 'Unit', 'normalized', 'Position', [0.06 0.25 0.8 0.2]);
    set(a2, 'Unit', 'normalized', 'Position', [0.06 0.49 0.8 0.2]);
    set(a3, 'Unit', 'normalized', 'Position', [0.06 0.73 0.8 0.2]);
    
%initial---
t = 0e-9;

[X Y] = meshgrid(16*geometry.dx:16*geometry.dx:2040*geometry.dx, 16*geometry.dy:16*geometry.dy:248*geometry.dy);
[X1 Y1] = meshgrid(8*geometry.dx:8*geometry.dx:2040*geometry.dx, 8*geometry.dy:8*geometry.dy:248*geometry.dy);

        load('-mat',strcat(path,'rho_electrons_',num2str(t),'.dat'),var) 
        rho_el = rho_sp;
%         load('-mat',strcat(path,'rho_ions_',num2str(t),'.dat'),var) 
%         rho_ions = rho_sp;     
        load('-mat',strcat(path,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;

        
%         rho_el = rho_el + rho_back*1.6e-19;

        
        rho = rho_el + rho_beam + rho_back;
        
        [ex ey fi] = field_2(rho, 1, 0, 0*ones(1,geometry.ngy), 0*ones(1,geometry.ngy)...
                      , 0*ones(1,geometry.ngx), 0*ones(1,geometry.ngx) ,geometry.dx, geometry.dy);
                  
        e_abs = (ex.^2 + ey.^2).^0.5;              
         
        q_ex = ex(16:16:248,16:16:2040);
        q_ey = ey(16:16:248,16:16:2040);
        c_fi = fi(8:8:248,8:8:2040);

        

%         imshow(rho,clim)
%         rho = imfilter(rho,h);

        im1 = image(rho_el,'Parent',a1, 'CDataMapping', 'scaled');

        
        set(a1, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [' 0.0'; '0.15'; ' 0.3'; '0.45'; ' 0.6'], ...
                'yTick', [128 256], 'yTickLabel', ['0.05'; ' 0.1'], 'Clim', clim1, 'YDir', 'normal','FontSize',14); 
        im2 = image(rho_ions,'Parent',a2, 'CDataMapping', 'scaled');
        set(a2, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim2, 'YDir', 'normal','FontSize',14);
%         im3 = image(rho_beam,'Parent',a3, 'CDataMapping', 'scaled');
%         set(a3, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim3, 'YDir', 'normal');

%        cont = contour(X1, Y1, c_fi, 'Parent',a3);
%        hold on;
%        quiv = quiver(a3, X,Y,q_ex,q_ey);
%        hold off;
       
       set(a3, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim3, 'YDir', 'normal','FontSize',14);
       
  

        
        %----------------------------------
    set(get(a1,'XLabel'),'String','x, m', 'FontSize', 16.0);
    set(get(a1,'YLabel'),'String','y, m', 'FontSize', 16.0);

    text_n0_0 = text(0,-0.5, 'n_i_,_e(x=0,t=0) = 2.4\cdot10^1^4 m^-^3', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
    text_n0_L = text( 0,-0.8, 'n_i_,_e(x=L,t=0) = 4.0\cdot10^1^4 m^-^3', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
    text_n0_L = text( 0,-1.1, 'T_e = 2.0 eV', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
    
    text_Te = text(0.3, -0.5, '\phi(x=0) = \phi(x=Lx) = 0; (Dirichlet bc)', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
    text_Ti = text(0.3, -0.8, '\phi(y=0) = \phi(y = Ly); (periodic bc)', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
    text_Ti = text(0.3, -1.1, 'T_i = 0.2 eV', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
    
    text_f = text(0.7, -0.5, 'f_m_o_d = 1.6\cdot10^8 Hz', 'FontSize', 14, 'Parent', a1,'Unit','normalized'); 
    text_n_b = text(0.7, -0.8, 'n_b_e_a_m = 1\cdot10^1^3 m^-^3', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
    text_v_b = text(0.7, -1.1, 'v_b_e_a_m = 50\cdotv_T', 'FontSize', 14, 'Parent', a1,'Unit','normalized');

    text_Nbe = text(0.95, -0.5, 'N_b_i_g_,_e_l = 1.2\cdot10^6', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
    text_Nbi = text(0.95, -0.8, 'N_b_i_g_,_i_o_n = 4.0\cdot10^5', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
    text_Ngr = text(0.95, -1.1, 'N_g_r_i_d = 2048\times256', 'FontSize', 14, 'Parent', a1,'Unit','normalized');
    
    text_title = text(0.4, 3.65, strcat('Beam-plasma interaction,  t =','  ',num2str(t/Tmod,'%2.1f'),'  T_m_o_d'), ...
    'FontSize', 16, 'Parent', a1,'Unit','normalized');
    
    text_n_fi = text(1.08, 0.9,'|E|, V/m' , 'Parent', a1,'Unit','normalized','FontSize', 16);
%     text_n_el_2 = text(1.1, 0.65, '10^1^4 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);
    
    text_n_ion_1 = text(1.1, 1.8, '|E|, averaged, V/m', 'Parent', a1,'Unit','normalized','FontSize', 16);
%     text_n_ion_2 = text(1.1, 1.85, '10^1^4 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);
    
%     text_n_beam_1 = text(1.1, 3.4,'\phi, V' , 'Parent', a1,'Unit','normalized','FontSize', 16);
%     text_n_beam_2 = text(1.1, 3.05, '10^1^2 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);

 
    clrbr1 = colorbar('peer', a1,'Position', [0.88 0.25 0.03 0.2], 'Clim',clim1,'XTick',[], 'xTickLabel', [],...            %
       'YTick',Tick_position_1, 'yTickLabel', Tick_label_1 );
    clrbr2 = colorbar('peer',a2, 'Position', [0.88 0.49 0.03 0.2],'Clim',clim2,'XTick',[], 'xTickLabel', [],...
        'YTick',Tick_position_2, 'yTickLabel', Tick_label_2 );
%     clrbr3 = colorbar('peer',a3, 'Position', [0.88 0.73 0.03 0.2],'Clim',clim3,'XTick',[], 'xTickLabel', [],...
%         'YTick',[-50 -100], 'yTickLabel', [' -50';'-100'] );
    
colormap('gray')

%----------
N_averaged_step = 15;
j = 1;
e_averaged = zeros(256,2048);
max_el = [];
for k = 1:100

    tstart = (k-1)*25*dt + t_begin;
    tend = (k*25-1)*dt + t_begin;
    %D = zeros(700,1050,3,75);
    i = 1;

    for t = tstart:dt:tend
        load('-mat',strcat(path,'rho_electrons_',num2str(t),'.dat'),var) 
        rho_el = rho_sp;    
        load('-mat',strcat(path,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;
        rho = rho_sp;
           
        rho = rho_el + rho_beam + rho_back;
        
        [ex ey fi] = field_2(rho, 1, 0, 0*ones(1,geometry.ngy), 0*ones(1,geometry.ngy)...
                      , 0*ones(1,geometry.ngx), 0*ones(1,geometry.ngx) ,geometry.dx, geometry.dy);
                  
        e_abs = (ex.^2 + ey.^2).^0.5;              
         
        q_ex = ex(16:16:248,16:16:2040);
        q_ey = ey(16:16:248,16:16:2040);
        c_fi = fi(8:8:248,8:8:2040);
        
        set(im1, 'CData', e_abs);
        
        [valy posy] = max(e_abs);
        [val posx] = max(valy);
        x = posx;
        y = posy(x);
        max_el(end+1) = val;
        
        plot('Parent', a3, max_el);
        
        hold on;
        plot('Parent', a1, x, y, 'MarkerFaceColor','g', 'MarkerSize',10);
        hold off;
        
        %----averaging
%             Averaging_array(:,:,2:10) = Averaging_array(:,:,1:9);
%             Averaging_array(:,:,1) = e_abs;
%             e_averaged = sum(Averaging_array,3)/10;
        %-------------
%         e_averaged = e_averaged + e_abs;
%         if j > N_averaged_step
%             t_back = t - N_averaged_step*dt;
%             load('-mat',strcat(path,'rho_electrons_',num2str(t_back),'.dat'),var) 
%             rho_el = rho_sp;
%             load('-mat',strcat(path,'rho_ions_',num2str(t_back),'.dat'),var) 
%             rho_ions = rho_sp;     
%             load('-mat',strcat(path,'rho_light_electrons_',num2str(t_back),'.dat'),var) 
%             rho_beam = rho_sp;
%             rho = rho_sp;
%         
%             rho_el = rho_el/1.6e-19 + rho_back;
%             rho_ions = rho_ions/1.6e-19 + rho_back*(-1);
%             rho_beam = rho_beam/1.6e-19;
%         
%             rho = (rho_el + rho_ions + rho_beam)*1.6e-19;
%         
%             [ex_back ey_back fi] = field_2(rho, 1, 0, 0*ones(1,geometry.ngy), 0*ones(1,geometry.ngy)...
%                       , 0*ones(1,geometry.ngx), 0*ones(1,geometry.ngx) ,geometry.dx, geometry.dy);
%                   
%             e_abs_back = (ex_back.^2 + ey_back.^2).^0.5;              
%             
%             e_averaged = e_averaged - e_abs_back;
%             
%             set(im2, 'CData', e_averaged/N_averaged_step);
%         else
%             set(im2, 'CData', zeros(256,2048));
%         end
            set(im2, 'CData', zeros(256,2048));      
        
%         set(im2, 'CData', fi);
     
%         cont = contour(X1, Y1, c_fi, 'Parent',a3, 'LineWidth', 2);
%        hold on;
%        set(a3, 'Unit', 'pixels');
%        quiv = quiver(a3, X,Y,q_ex/800*8/2048*0.6,q_ey/800*8/256*0.1, 'Autoscale', 'off','LineWidth',1.2, 'MaxHeadSize', 1.7);
%        hold off;
% %        set(a3, 'Unit', 'normalized');
%        set(a3, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim3, 'YDir', 'normal','FontSize',14);
%         
%         set(a2,'xlim',[0 0.6])


set(text_title, 'String', strcat('Beam-plasma interaction,  t =','  ',num2str(t/Tmod,'%2.1f'),'  T_m_o_d'));


        figHandle = gcf;
        frame = getframe(figHandle);
        D(:,:,:,i)   = frame.cdata;
        D(:,:,:,i+1) = frame.cdata;
        D(:,:,:,i+2) = frame.cdata;
       
        i = i + 3;
        j = j + 1;
    end

    
    f2 = figure('Position', [20 20 100 100]);
    at = axes('Parent', f2);
    mov = immovie(D);
    movie2avi(mov, strcat(path4wr,'el_movie_',num2str(t/1e-7,'%3.2f'),'.avi'));
    clear D mov
    close(f2)
    pack



end
res = 1;
