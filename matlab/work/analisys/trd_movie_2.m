function res = trd_movie_2;

dt = 5e-10;

w = 1.0e9;
Tmod = 2*pi/w;

h = ones(10,10)/100;


path = 'D:\_results\mobile_1.6e14_10e12_01Ly_05.31_overcrit4\';
load('D:\_results\mobile_1.6e14_10e12_01Ly_05.31_overcrit4\mobile_4.8e14_10e12_01Ly_overcrit4.mat', 'param');

geometry = param.geometry;
bc = param.bc;
geometry = calc_grid_step(geometry,bc);
rho_back_struct = param.rho_back_struct;


t_begin = 2.0e-8;
var = 'rho_sp';


if rho_back_struct.enabled
    rho_back = load_rho_back(rho_back_struct, geometry, bc);
end



%-----determine visualization range-------


        clim2 = [-1.0e13 0];
        clim3 = [-4.0e14 4.0e14];
        clim1 = [-50 50];
   
  


%----------
h = ones(10)/100;
h2= ones(20)/4;

% 
%     tstart = (k-1)*25*dt + t_begin;
%     tend = (k*25-1)*dt + t_begin;
%     i = 1;
%     tstart = 6.54e-7;
%     tend = 6.78e-7;

%     e_av = zerros(256,2048);
%     N = length(tstart:dt:tend);
%     for t = tstart:dt:tend
        t = 2.0e-8;
        load('-mat',strcat(path,'rho_electrons_',num2str(t),'.dat'),var) 
        rho_el = rho_sp;
%         load('-mat',strcat(path,'rho_ions_',num2str(t),'.dat'),var) 
%         rho_ions = rho_sp;     
        load('-mat',strcat(path,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;
        rho = rho_sp;
        rho = rho_el+rho_back*1.6e-19+rho_beam;
        
        rho_el = rho_el/1.6e-19 + rho_back;
        rho_el = rho_el*(-1);
%         rho_ions = rho_ions/1.6e-19 + rho_back*(-1);
        rho_beam = rho_beam/1.6e-19;
        
        fi = field_3(rho,geometry,bc);
        [ex ey] = e_from_fi(fi,geometry,bc);
%         e_av = (ex.^2 + ey.^2).^0.5 + e_av;
%     end
%     e_av = e_av/N;
%        set(im1, 'CData', fi);
%        set(im2, 'CData', rho_beam);

%        p1 = plot(fi(128,:), 'Parent', a3);
%        p2 = plot(sum(rho_beam), 'Parent', a2);
      

  f = figure('Position',[20 100 1200 400]);
  a = axes('Parent', f, 'Color', 'none', 'XTick',[0 0.2 0.4 0.6],'XTickLabel',['0.0';'0.2';'0.4';'0.6'],...
                                         'YTick',[0 0.05 0.1], 'YTickLabel', [' 0.0';'0.05';' 0.1'],...
                                         'clim', [-100 100],...
                                         'ZLim', [-100 100]...
                                     );
  xlabel('x, m','FontSize',16)
  ylabel('y, m','FontSize',16)
  zlabel('\phi, V','FontSize',16)
%   hold on;
[x,y] = meshgrid(0:geometry.dx:geometry.x_size, 0:geometry.dy:geometry.y_size);

% fi = fi(1:2:end,:);
% rho_beam = rho_beam(1:2:end,:);
% x = x(1:2:end,:);
% y = y(1:2:end,:);
  s = surf(x,y,imfilter(fi,h), 'Parent', a, 'LineStyle', 'none');
%   AData = rho_beam/max(max(abs(rho_beam)))*(-1);
%   hold on;
%   s = surf(imfilter(fi,h), 'Parent', a, 'LineStyle', 'none', 'Cdata', rho_beam, 'CDataMapping', 'scaled',...
%       'FaceAlpha', 'flat', 'AlphaData', rho_beam, 'AlphaDataMapping','scaled');
%   hold off;

  hold on;
  surf(x,y,imfilter(fi,h),'Parent',a,'LineStyle','none',...
     'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',imfilter(rho_beam*(-1),h2),...
    'FaceColor','blue');
axis tight

zlim([-100 100])
view([10 70])    
colormap('copper')



        
%     end


res = rho_beam;

end

