function res = draw_plot(input);


filter = ones(15)/225;

% input = imfilter(input,filter);

clim1 = [min(min(input))/1.5 max(max(input))/1.1];

% clim1 = [0 2];

f = figure;
set(f, 'Position', [50 60 1050 700], 'Color', 'white');    
 
a1 = axes('Parent', f);

set(a1, 'Unit', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);


        
im1 = image(input,'Parent',a1, 'CDataMapping', 'scaled');

set(a1, 'xTick', [1 512 1024 1536 2048], 'xTickLabel', [' 0.0'; '0.15'; ' 0.3'; '0.45'; ' 0.6'],...
    'yTick', [1 500 1000 1500], 'yTickLabel', ['0.0'; '0.5'; '1.0';'1.5'], 'Clim', clim1, 'YDir', 'normal'); 


res = 1;
