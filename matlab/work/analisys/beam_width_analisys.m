function [beam_width_3d beam_density_3d] = beam_width_analisys


% n = 2e12;
% path = 'G:\_results\mobile_24-40e14_20e12_01Ly_04.20\';
% load('G:\_results\mobile_24-40e14_20e12_01Ly_04.20\mobile_24-40e14_20e12_01Ly.mat', 'param');

% nb = 1e12; L_inhom = 1*Lx;
% path = 'G:\_results\mobile_16-48e14_05e11_01Ly_04.02\';
% load('G:\_results\mobile_16-48e14_05e11_01Ly_04.02\mobile_16-48e14_05e11_01Ly.mat', 'param');

% nb = 1e12; L_inhom = 1.5*Lx;
% path = 'G:\_results\mobile_20-44e14_05e11_01Ly_04.05\';
% load('G:\_results\mobile_20-44e14_05e11_01Ly_04.05\mobile_20-44e14_05e11_01Ly.mat', 'param');

% path = 'G:\_results\mobile_16e14_10e12_01Ly_08.17_subcrit\';
% load('G:\_results\mobile_16e14_10e12_01Ly_08.17_subcrit\mobile_16e14_10e12_01Ly_subcrit.mat', 'param');

% path = 'G:\_results\mobile_48e14_10e12_01Ly_08.20_overcrit\';
% load('G:\_results\mobile_48e14_10e12_01Ly_08.20_overcrit\mobile_48e14_10e12_01Ly_overcrit.mat', 'param');

path = 'G:\_results\bd_48e14_10e12_01Ly_periodic_10_09_supercrit\';
load('G:\_results\bd_48e14_10e12_01Ly_periodic_10_09_supercrit\bd_48e14_10e12_01Ly_periodic_supercrit.mat', 'param');

% nb = 1e12; L_inhom = 2.0*Lx;
% path = 'G:\_results\mobile_24-40e14_05e11_01Ly_04.08\';
% load('G:\_results\mobile_24-40e14_05e11_01Ly_04.08\mobile_24-40e14_05e11_01Ly.mat', 'param');

geometry = param.geometry;
bc = param.bc;
rho_back_struct = param.rho_back_struct;

geometry = calc_grid_step(geometry, bc);

init_beam_width = param.beam_struct(1).inject_param.beam_width;

dt = 5e-10;

w = 1.0e9;
Tmod = 2*pi/w;


var = 'rho_sp';

filter_length = 50;
filter = ones(1,filter_length);


t_begin = 0;
%----------
t_start = 0.5e-7;
t_end = 3.2e-7;


const.eq = 1.6e-19;
if rho_back_struct.enabled
    rho_back = load_rho_back(rho_back_struct, geometry, bc)*const.eq;
end


k = 1;    
beam_width_3d = zeros(length(t_start:dt:t_end), geometry.ngx - 2*filter_length);
beam_density_3d = zeros(length(t_start:dt:t_end), geometry.ngx - 2*filter_length);


for t = t_start:dt:t_end
    
            try,
        load('-mat',strcat(path,'rho_light_electrons_',num2str(t),'.dat'),var) 
        rho_beam = rho_sp;
        rho = rho_sp;
        catch,
            
            ttt = t+dt/10;
        load('-mat',strcat(path,'rho_light_electrons_',num2str(ttt),'.dat'),var) 
        rho_beam = rho_sp;
        rho = rho_sp;
        end
    

        n_beam_1d = sum(rho_beam)/geometry.ngy/const.eq/init_beam_width;


        for i = 1:2048
            vector1 = find(rho_beam(:,i)' ~= 0);

            if length(vector1) > 0
                L1(i) = vector1(1);
            else 
                L1(i) = 128;
            end
                
            vector2 = find(fliplr(rho_beam(:,i)') ~= 0);
            if length(vector2) > 0
                L2(i) = vector2(1);
            else
                L2(i) = 128;
            end
        end
        width = max(abs(L2-128),abs(L1 - 128));

       
n_beam_1d = imfilter(n_beam_1d, filter)/filter_length;
width = imfilter(width, filter)/filter_length*2*geometry.y_size/geometry.ngy;

n_beam_1d = n_beam_1d(filter_length+1:geometry.ngx-filter_length);
width = width(filter_length+1:geometry.ngx-filter_length);


%interpolation
%x_start  the start point of interpolation interval
%x_end - the end point of interpolation interval
%x_low_density_start - the start point of low density region
%x_low_density_end - the end point of low density region


max_density_value = max(abs(n_beam_1d))
threshhold = 0.05;

% rare_n_beam_1d = n_beam_1d(5:50:end);


% need_correction_vector = find(n_beam_1d > (-1)*max_density_value*threshhold);
%  
% no_need_correction_vector = find(n_beam_1d < (-1)*max_density_value*threshhold);
% 
% width1 = interp1(no_need_correction_vector, width(no_need_correction_vector), 1:length(width), 'spline');

need_correction_vector = find(n_beam_1d > (-1)*max_density_value*threshhold);
 
no_need_correction_vector = find(n_beam_1d < (-1)*max_density_value*threshhold);


width1 = interp1(no_need_correction_vector(5:1:end), width(no_need_correction_vector(5:1:end)), 1:length(width), 'spline');


% width1 = imfilter(width1(10:end), filter)/filter_length*2*geometry.y_size/geometry.ngy;
% 
% 
% width1 = width1(filter_length+1:geometry.ngx-filter_length);

% figure, plot(width1)
% adf

j = 1;
x_start = 0;
x_end = 0;

% for i = need_correction_vector
% 
%     if (i < x_end)&(i > x_start)
%         if x_start == 0
%             width(i) = width(x_end);
%         elseif x_end == geometry.ngx - 2*filter_length +1
%             width(i) = width(x_start);
%         else
%             width(i) = width(x_start)*(x_end - i)/(x_end - x_start) + width(x_end)*(i - x_start)/(x_end - x_start);
%         end
%     else
%         x_start = need_correction_vector(j) - 1;
% 
%         while true
%             if need_correction_vector(j+1) - need_correction_vector(j) < 10
%                 if j < length(need_correction_vector) - 1
%                     j = j + 1;
%                     continue;
%                 else
%                     x_end = need_correction_vector(end) + 1;
%                     break;
%                 end
%             else
%                 x_end = need_correction_vector(j) + 1;
%                 if j < length(need_correction_vector) - 1
%                     j = j + 1;
%                 end
%                 break;
%             end        
%         end
% %         i = x_start + 1;
%         
%     end
%     
% end

beam_width_3d(k, :) = width1;
beam_density_3d(k, :) = n_beam_1d;
k = k + 1


% figure        
% subplot(2,1,1), plot(n_beam_1d)
% subplot(2,1,2), plot(width)
% afdd


end




