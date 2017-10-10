function res = compare_two_simulations 


dt = 5e-10;
w = 1.0e9;
Tmod = 2*pi/w;

load('d:\_results\mobile_20-44e14_10e12_01Ly_04.23\mobile_20-44e14_10e12_01Ly.mat');
geometry = param.geometry;
bc = param.bc;

geometry = calc_grid_step(geometry, bc);


rho_back_struct(1) = struct('name', {'ions'},  ...
                            'charge', {1},     ...
                            'profile_fun',  {struct('handle',{@linear_fun}, 'param',{[2.4e14 4.0e14]})},...
                            'enabled', {1});



h = ones(5,5)/25;

path_1 = 'g:\_results\mobile_20-44e14_05e11_01Ly_04.05\';
path_2 = 'd:\_results\mobile_20-44e14_10e12_01Ly_04.23\';
path4wr = 'd:\_results\movies\';
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

        clim1 = [0 12e3];
        clim4 = [0 12e3];
        clim3 = [-5e12 0];
        clim2 = [-5e12 0];
        time_lim = [0 150];
        
        Tick_position_1 = [4000 8000];
        Tick_position_4 = [4000 8000];
        Tick_position_3 = [-2E-12 -4E-12];
        Tick_position_2 = [-2E-12 -4E-12];
        
        Tick_label_1 = ['4000'; '8000'];
        Tick_label_4 = ['4000'; '8000'];
        Tick_label_3 = ['-2e-12'; '-4e-12'];
        Tick_label_2 = ['-2e-12'; '-4e-12'];
  
f = figure;
set(f, 'Position', [50 60 1050 700], 'Color', 'white');    
 
a1 = axes('Parent', f);
a2 = axes('Parent', f);
a3 = axes('Parent', f);
a4 = axes('Parent', f);
    set(a1, 'Unit', 'normalized', 'Position', [0.06 0.25 0.8 0.2], 'NextPlot','add');
    set(a2, 'Unit', 'normalized', 'Position', [0.06 0.49 0.8 0.2]);
    set(a3, 'Unit', 'normalized', 'Position', [0.06 0.73 0.8 0.2]);
    set(a4, 'Unit', 'normalized', 'Position', [0.06 0.02 0.8 0.2]);
    
%initial---
       

%         im2 = image(zeros(256,2048),'Parent',a2, 'CDataMapping',
%         'scaled');
    zero = zeros(256,2048);
    im1 = image(zero,'Parent',a1, 'CDataMapping', 'scaled');    
    im2 = image(zero,'Parent',a2, 'CDataMapping', 'scaled');  
    im3 = image(zero,'Parent',a3, 'CDataMapping', 'scaled');  
    im4 = image(zero,'Parent',a4, 'CDataMapping', 'scaled');      
    
    set(a1, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim1, 'YDir', 'normal'); 
    set(a2, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim2, 'YDir', 'normal'); 
    set(a3, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim3, 'YDir', 'normal'); 
    set(a4, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [], 'yTick', [128 256], 'yTickLabel', [], 'Clim', clim4, 'YDir', 'normal'); 
        %----------------------------------

    
    text_title = text(0.4, 3.65, strcat('Beam-plasma interaction,  t =','  ',num2str(t_begin/Tmod,'%2.1f'),'  T_m_o_d'), ...
    'FontSize', 16, 'Parent', a1,'Unit','normalized');
    
    text_axes_1 = text(1.08, 0.9,'|E|1, mob' , 'Parent', a1,'Unit','normalized','FontSize', 16);
    text_axes_4 = text(1.08, -0.9,'|E|2, imm' , 'Parent', a1,'Unit','normalized','FontSize', 16);

    
    text_axes_21 = text(1.065, 2.15, 'Beam2, imm', 'Parent', a1,'Unit','normalized','FontSize', 16);
%     text_axes_22 = text(1.065, 1.95, 'averaged', 'Parent', a1,'Unit','normalized','FontSize', 16);
    
    text_axes_3 = text(1.08, 3.4, 'Beam1, mob', 'Parent', a1,'Unit','normalized','FontSize', 16);
%     text_n_ion_2 = text(1.1, 1.85, '10^1^4 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);
    
%     text_n_beam_1 = text(1.1, 3.4,'\phi, V' , 'Parent', a1,'Unit','normalized','FontSize', 16);
%     text_n_beam_2 = text(1.1, 3.05, '10^1^2 m^-^3', 'Parent', a1,'Unit','normalized','FontSize', 12);

 
    clrbr1 = colorbar('peer', a1,'Position', [0.88 0.25 0.03 0.2], 'Clim',clim1,'XTick',[], 'xTickLabel', [],...            %
       'YTick',Tick_position_1, 'yTickLabel', Tick_label_1 );
    clrbr2 = colorbar('peer',a2, 'Position', [0.88 0.49 0.03 0.2],'Clim',clim2,'XTick',[], 'xTickLabel', [],...
        'YTick',Tick_position_2, 'yTickLabel', Tick_label_2 );
    clrbr3 = colorbar('peer',a3, 'Position', [0.88 0.73 0.03 0.2],'Clim',clim3,'XTick',[], 'xTickLabel', [],...
        'YTick',Tick_position_3, 'yTickLabel', Tick_label_3 );
    clrbr4 = colorbar('peer',a4, 'Position', [0.88 0.02 0.03 0.2],'Clim',clim4,'XTick',[], 'xTickLabel', [],...
        'YTick',Tick_position_4, 'yTickLabel', Tick_label_4 );

colormap('gray')

%----------

for k = 1:100

    tstart = (k-1)*25*dt + t_begin;
    tend = (k*25-1)*dt + t_begin;
    i = 1;

    for t = tstart:dt:tend
        %----------first-----simulation
        load('-mat',strcat(path_1,'rho_electrons_',num2str(t),'.dat'),var) 
        rho_el = rho_sp;    
        load('-mat',strcat(path_1,'rho_ions_',num2str(t),'.dat'),var) 
        rho_ions = rho_sp;   
        load('-mat',strcat(path_1,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;
        rho = rho_sp;
        
        rho = rho_el + rho_beam + rho_ions;
        
        fi = field_3(rho, geometry, bc);
        [ex ey] = e_from_fi(fi, geometry, bc);          
        e_abs = (ex.^2 + ey.^2).^0.5;              
                
        rho_beam = rho_beam/1.6e-19;
        set(im3, 'CData', rho_beam);
        
        set(im1, 'CData', e_abs);
  
        [valy posy] = max(e_abs);
        [val posx] = max(valy);
        x1 = posx;
        y1 = posy(x1);

         
                
        %---------second---------beam
        
        load('-mat',strcat(path_2,'rho_electrons_',num2str(t),'.dat'),var) 
        rho_el = rho_sp;      
        load('-mat',strcat(path_1,'rho_ions_',num2str(t),'.dat'),var) 
        rho_ions = rho_sp; 
        load('-mat',strcat(path_2,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;
        rho = rho_sp;
        
        rho = rho_el + rho_beam + rho_ions;
        

        
        fi = field_3(rho, geometry, bc);
        [ex ey] = e_from_fi(fi, geometry, bc);          
        e_abs = (ex.^2 + ey.^2).^0.5;              
                
        rho_beam = rho_beam/1.6e-19;
        set(im2, 'CData', rho_beam);
        
        set(im4, 'CData', e_abs);

        [valy posy] = max(e_abs);
        [val posx] = max(valy);
        x2 = posx;
        y2 = posy(x2);
        
        
%        %-------------------------------        
        p1 = plot(x1, y1,'-.or', 'MarkerFaceColor','w','MarkerEdgeColor','k', 'MarkerSize',10, 'Parent', a1);
                set(a1, 'ylim',[1 256], 'xlim',[1 2048])
%         p2 = plot(x2, y2,'-.or', 'MarkerFaceColor','w','MarkerEdgeColor','k', 'MarkerSize',10, 'Parent', a4);
%                 set(a4, 'ylim',[1 256], 'xlim',[1 2048])




set(text_title, 'String', strcat('Beam-plasma interaction,  t =','  ',num2str(t/Tmod,'%2.1f'),'  T_m_o_d'));


        figHandle = gcf;
        frame = getframe(figHandle);
        D(:,:,:,i)   = frame.cdata;
        D(:,:,:,i+1) = frame.cdata;
        D(:,:,:,i+2) = frame.cdata;
        whos 'D'
        i = i + 3;
        delete(p1)
%         delete(p2)
    end

    
    f2 = figure('Position', [20 20 100 100]);
    at = axes('Parent', f2);
    mov = immovie(D);
    movie2avi(mov, strcat(path4wr,'el_movie_',num2str(t/1e-7,'%3.2f'),'.avi'));
    %clear mov
    close(f2)
    %pack



end
res = 1;
