function var = ions_movie_create;

w = 1.0e9;
Tmod = 2*pi/w;

filter = ones(9,12)/9/12;


path = 'G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\';
path4wr = 'G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\_movie\';
load('G:\_results\ion_motion_16-48_10e12_01Ly_periodic_05_12\ion_motion_16-48_10e12_01Ly_periodic.mat', 'param');

geometry = param.geometry;
bc = param.bc;

geometry = calc_grid_step(geometry, bc);
rho_back_struct = param.rho_back_struct;

t_begin = 0.0e-7;
dt = 5e-9;
var = 'rho_sp';

% if rho_back_struct.enabled
    rho_back = load_rho_back(rho_back_struct, geometry, bc)*1.6e-19;

% end0

%-----determine visualization range-------



%         clim1 = [-3e14 3e14];
%         clim2 = [-3e14 3e14];
%         clim3 = [-0.6e13 0.6e13];
%         
%         tick1 = [-2.0E14 0 2.0E14];
%         tick2 = [-2.0E14 0 2.0E14];
%         tick3 = [-4E12 0 4E12];
%         
%         tick_label_1 = ['-2'; ' 0'; ' 2'];
%         tick_label_2 = ['-2'; ' 0'; ' 2'];
%         tick_label_3 = ['-4';' 0'; ' 4'];        
        
        clim1 = [0 6e3];
%         clim2 = [-300 300];
%         clim2 = [-3e14 3e14];        
        clim2 = [1e14 5e14];        
        clim3 = [-6e12 6e12];
        
        tick1 = [2e3 4e3 6e3];
%         tick2 = [-200 0 200];
%         tick2 = [-2.0E14 0 2.0E14];        
        tick2 = [-2.0E14 0 2.0E14 3.2E14];        
        tick3 = [-4E12 0 4E12];
        
        tick_label_1 = ['2'; '4'; '6'];
%         tick_label_2 = ['-200'; '   0'; ' 200'];
%         tick_label_2 = ['-2'; ' 0'; ' 2'];
        tick_label_2 = [' -2'; '  0'; '  2'; '3.2'];
        tick_label_3 = ['-4';' 0'; ' 4'];        
   
  
f = figure;

set(f, 'Position', [50 60 1100 200], 'Color', 'white');  

a1 = axes('Parent', f);
set(a1, 'Unit', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);

z = zeros(geometry.ngy, geometry.ngx);
im1 = image(z,'Parent',a1, 'CDataMapping', 'scaled');

grid_points = [1 512 1024 1536 2048];
labels = [' 0.0'; '0.15'; ' 0.3'; '0.45'; ' 0.6'];


grid_points = [1    512     582     652     722     792     862     932      1002        1094 1536 2048];
labels = [' 0.0'; '0.15'; '0.17'; '0.19'; '0.21'; '0.23'; '0.25'; '0.27'; '0.29';  '0.45'; ' 0.6']


set(a1, 'xTick', grid_points, 'xTickLabel', labels,...
    'yTick', [128 256], 'yTickLabel', ['0.05'; ' 0.1'], 'Clim', [-3e14 3e14], 'YDir', 'normal'); 

text_1 = text( 1.01,0.5, 'T=', 'Parent', a1,'Unit','normalized','FontSize', 16);


for k = 1:150

    tstart = (k-1)*50*dt + t_begin;
    tend = (k*50-1)*dt + t_begin;
    i = 1;

    for t = tstart:dt:tend
        
        try,
        load('-mat',strcat(path,'rho_electrons_',num2str(t),'.dat'),'rho_sp') 
        rho_el = rho_sp;

        load('-mat',strcat(path,'rho_ions_',num2str(t),'.dat'),var) 
        rho_ions = rho_sp;     
        load('-mat',strcat(path,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;
        catch,
            
            ttt = t+dt/20/5;
                load('-mat',strcat(path,'rho_electrons_',num2str(ttt),'.dat'),'rho_sp') 
        rho_el = rho_sp;

        load('-mat',strcat(path,'rho_ions_',num2str(ttt),'.dat'),var) 
        rho_ions = rho_sp;     
        load('-mat',strcat(path,'rho_light_electrons_',num2str(ttt),'.dat'),var) 
        rho_beam = rho_sp;
        end
            

%         rho = rho_sp;
        rho = rho_el+rho_ions+rho_beam;
%         rho = rho_el + rho_beam + rho_back;
        
        
%         rho_el = (rho_el + rho_back)/1.6e-19;
%         rho_el = rho_el*(-1);
%         rho_ions = rho_ions/1.6e-19 + rho_back*(-1);
        rho_beam = rho_beam/1.6e-19;
        
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

        
        
 
%         set(im1, 'CData', sum(ex_av,3)/N);        
%         set(im2, 'CData', sum(ey_av,3)/N);
%         set(im3, 'CData', ((sum(ex_av,3)).^2 + (sum(ex_av,3)).^2).^0.5/N);
        
       rho_el_delta = (rho_el + rho_back)/1.6e-19*(-1);
       rho_ion_delta = (rho_ions - rho_back)/1.6e-19;
       
% 
        set(im1, 'CData', imfilter(rho_ion_delta,filter));
%         imshow(imfilter(rho_ion_delta,filter),[-3e14 3e14]);

set(text_1, 'String',strcat('f_p\cdott =  ','  ',num2str(t/Tmod,'%2.1f')));
colormap('jet')

        figHandle = gcf;
        frame = getframe(figHandle);
        D(:,:,:,i)   = frame.cdata;

        i = i + 1;
   drawnow;
        
    end
    f2 = figure('Position', [20 20 100 100]);
    at = axes('Parent', f2);
    mov = immovie(D);
    movie2avi(mov, strcat(path4wr,'rho_ions_movie_2_',num2str(t/1e-7,'%3.2f'),'.avi'));
    clear D mov
    close(f2)
    pack



end
res = 1;
