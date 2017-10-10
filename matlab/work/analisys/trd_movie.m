function res = trd_movie;

w = 1.0e9;
Tmod = 2*pi/w;

path = 'D:\_results\mobile_3.2e14_10e12_01Ly_06_05_critical2\';
path4wr = 'D:\_results\mobile_3.2e14_10e12_01Ly_06_05_critical2\_movie\';
load('D:\_results\mobile_3.2e14_10e12_01Ly_06_05_critical2\mobile_3.2e14_10e12_01Ly_critical.mat', 'param');

geometry = param.geometry;
bc = param.bc;

geometry = calc_grid_step(geometry, bc);
rho_back_struct = param.rho_back_struct;


dt = 5e-10;
var = 'rho_sp';

h = ones(10)/100;
h2 = ones(20)/4;

% if rho_back_struct.enabled
    rho_back = load_rho_back(rho_back_struct, geometry, bc)*1.6e-19;
% end

%-----determine visualization range-------

clim1 = [-150 150];

tick_x = [0.0 0.2 0.4 0.6];
tick_y = [0.0 0.05 0.1];
        
tick_label_x = ['0.0'; '0.2'; '0.4'; '0.6'];
tick_label_y = [' 0.0'; '0.05'; ' 0.1'];

fig_width = 1200;
fig_height = 400;

  f = figure('Position',[20 100 fig_width fig_height]);
  a = axes('Parent', f, 'Color', 'none', 'XTick',[0 0.2 0.4 0.6],'XTickLabel',['0.0';'0.2';'0.4';'0.6'],...
                                         'YTick',[0 0.05 0.1], 'YTickLabel', [' 0.0';'0.05';' 0.1'],...
                                         'clim', [-100 100],...
                                         'ZLim', [-100 100],...
                                         'FontSize',16 ...
                                     );
  xlabel('x, m','FontSize',16)
  ylabel('y, m','FontSize',16)
  zlabel('-\phi, V','FontSize',16)

t_begin = 0.0e-7;  

%initial-------------------------
[x,y] = meshgrid(0:geometry.dx:geometry.x_size, 0:geometry.dy:geometry.y_size);
  s1 = surf(x,y,zeros(256,2048), 'Parent', a, 'LineStyle', 'none');

  hold on;
  s2 = surf(x,y,zeros(256,2048),'Parent',a,'LineStyle','none',...
     'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',zeros(256,2048),...
    'FaceColor','blue');
%  set(a,'CLim',[-1e13 -1e11]);
axis tight
hold off;


zlim(clim1)
view([10 80]) 
colormap('copper')


% cl = camlight('left')
%--------------------------------
N = 5;
% D = ones(fig_height, fig_width, 3, N*3);
  
for k = 1:200

    tstart = (k-1)*N*dt + t_begin;
    tend = (k*N-1)*dt + t_begin;
    i = 1;

    for t = tstart:dt:tend
        load('-mat',strcat(path,'rho_electrons_',num2str(t),'.dat'),var) 
        rho_el = rho_sp;
        rho_el_1 = rho_sp;
%         load('-mat',strcat(path,'rho_ions_',num2str(t),'.dat'),var) 
%         rho_ions = rho_sp;     
        load('-mat',strcat(path,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;

%         rho = rho_el+rho_ions+rho_beam;
        rho = rho_el + rho_beam + rho_back;
              
        fi = field_3(rho,geometry,bc);
%         [ex ey] = e_from_fi(fi,geometry,bc);
%         e_abs = (ex.^2 + ey.^2).^0.5;
        rho_beam = rho_beam/1.6e-19;
        title(strcat('Beam-plasma interaction,  t =','  ',num2str(t/Tmod,'%2.1f'),'  T_m_o_d', ', critical plasma'), 'FontSize', 16);




%   s1 = surf(x,y,imfilter(fi,h), 'Parent', a, 'LineStyle', 'none');
set(s1, 'ZData', imfilter(fi,h)*(-1));

  hold on;
%   s2 = surf(x,y,imfilter(fi,h),'Parent',a,'LineStyle','none',...
%      'FaceAlpha','flat',...
%     'AlphaDataMapping','scaled',...
%     'AlphaData',imfilter(rho_beam*(-1),h2),...
%     'FaceColor','blue');
set(s2, 'ZData', imfilter(fi,h)*(-1), 'AlphaData', (rho_beam*(-1)).^0.5);
axis tight
hold off;
set(a, 'CLim', clim1, 'ZLim', clim1);


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
    movie2avi(mov, strcat(path4wr,'3D_movie_',num2str(t/1e-7,'%3.2f'),'.avi'));

    close(f2)
%     pack



end
res = 1;
