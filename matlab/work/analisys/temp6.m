function res = temp6

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

% rho_back_struct(1) = struct('name', {'ions'},  ...
%                             'charge', {1},     ...
%                             'profile_fun',  {struct('handle',{@linear_fun}, 'param',{[1.6e14 4.8e14]})},...
%                             'enabled', {1});

rho_back_struct(1) = struct('name', {'ions'},  ...
                            'charge', {1},     ...
                            'profile_fun',  {struct('handle',{@linear_fun}, 'param',{[2.4e14 4.0e14]})},...
                            'enabled', {0});

% xSize = 2e-0;
% ySize = 0.5e-0;
% 
% ngx = 512;
% ngy = 128;

% h = ones(5,5)/25;

path = 'd:\_results\moveless_24-40e14_10e12_01Ly_04.15\';

% path4wr = 'd:\temp\11.15\movie\';
% path = 'c:\movie';
var = 'rho_sp';
% name_sp = 'rho_beam_';

% f = figure;

t_begin = 0.00e-7;

%----------
t_start = 0.0e-7 + t_begin;
t_end = 7.0e-7 + t_begin;

% load('d:\_results\04.10\planar_waves_2L_05e11_wide.mat', 'param')
load('d:\_results\moveless_24-40e14_10e12_01Ly_04.15\moveless_24-40e14_10e12_01Ly.mat', 'param')
geometry = param.geometry;
if strcmp(bc.x_type, 'periodic')
    geometry.dx = geometry.x_size/geometry.ngx;
else
    geometry.dx = geometry.x_size/(geometry.ngx - 1);
end

if strcmp(bc.y_type, 'periodic')
    geometry.dy = geometry.y_size/geometry.ngy;
else
    geometry.dy = geometry.y_size/(geometry.ngy - 1);
end
bc = param.bc;
const.eq = 1.6e-19;
rho_back = load_rho_back(rho_back_struct, geometry, bc)*const.eq;

% 
% for k = 1:100
% 
%     tstart = (k-1)*25*dt + t_begin;
%     tend = (k*25-1)*dt + t_begin;
    i = 1;
    e_max_value_vector = zeros(size(t_start:dt:t_end));

    for t = t_start:dt:t_end
       load('-mat',strcat(path,'rho_electrons_',num2str(t),'.dat'),var) 
        rho_el = rho_sp;
%         load('-mat',strcat(path,'rho_ions_',num2str(t),'.dat'),var) 
%         rho_ions = rho_sp;     
        load('-mat',strcat(path,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;
       
%         rho = rho_beam + rho_el;
%         rho_beam_filtered = imfilter(rho_beam,ones(5)/25);
        
        rho = (rho_el + rho_back + rho_beam);
        


        
%         rho_beam = ones(256,2048)*1e11*1.6e-19;
%         rho = rho_ions;
%         [ex ey fi] = field_2(rho, 1, 0, 0*ones(1,geometry.ngy), 0*ones(1,geometry.ngy)...
%                       , 0*ones(1,geometry.ngx), 0*ones(1,geometry.ngx) ,geometry.dx, geometry.dy);
 
        fi = field_3(rho, geometry, bc);
        [ex ey] = e_from_fi(fi, geometry, bc);
        
        
        e_abs = (ex.^2 + ey.^2).^0.5;
        
                
%         figure, imshow(rho_el+rho_back,[])
%         figure, imshow(rho_ions-rho_back,[])
%         figure, plot(sum(rho_ions(96:160,:)))
        
%         figure, imshow(ex,[])

        
%         figure, surf(ex)


%         res = sum(ex)/256;
%         return;
        
%         [valy posy] = max(e_abs);
%         [val posx] = max(valy);
%         x = posx
%         y = posy(x)
        
        e_max_value_vector(i) = max(max(e_abs(:,700:1400)));
        plot(e_max_value_vector);
        drawnow;
         i = i+1

        
    end


% end
res = e_max_value_vector;
