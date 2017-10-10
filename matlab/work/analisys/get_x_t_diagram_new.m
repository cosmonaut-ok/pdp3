function res = get_x_t_diagram;

filter = ones(3,2)/6;

path = 'H:\PDP2\11.09\';
load('H:\PDP2\11.09\config_real.mat', 'param');

t_0 = 0;
dt = 0.25e-9;
t_max = 1.9505e-008;

x_t_diagram = zeros(length(t_0:dt:t_max), param.geometry.ngx);


i=1;
for t=t_0:dt:t_max
        try,
            load('-mat', strcat(path, 'rho_plasma_electrons_', num2str(t), '.dat'), 'rho_sp') 
        catch,
            t_new = t + dt/50;
            load('-mat', strcat(path, 'rho_plasma_electrons_', num2str(t_new), '.dat'), 'rho_sp')
        end
        rho_sp_filter = imfilter(rho_sp, filter);

        x_t_diagram(i,:) = rho_sp_filter(128,:) + 1e17*1.6e-19;

        i = i + 1;
end

res = x_t_diagram;

f = figure;
a = axes('Parent', f);
im1 = image(flipud(res),'Parent',a, 'CDataMapping', 'scaled');
colormap('gray');