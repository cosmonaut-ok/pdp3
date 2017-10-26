% rho_movie_create_light3('~/pdp3_modeling/model5/pdp3_result/', '~/pdp3_modeling/model5/movie/',100,254,2046,[0 1],[0 1],[-1e-7 0])

function var = rho_movie_create_light3(data_path,
                                       video_path,
                                       file_delta=100,
                                       size_1=254,
                                       size_3=2046,
                                       clim1=[0 1],
                                       clim2=[0 1],
                                       clim3=[-1e-7 0]);

                                % ,100,254,2046,[0 1],[0 1],[-1e-7 0])

data_file1 = strcat(data_path, '/', 'e1');
data_file2 = strcat(data_path, '/', 'e3');
data_file3 = strcat(data_path, '/', 'rho_beam');

% w = 1.0e9;
% Tmod = 2*pi/w;

filter = ones(3,12)/3/12;

path4wr = video_path;
name4wr = 'h2';

t_begin = 1.6e-6;
dt = 1e-9;
var = 'rho_sp';

f = figure;
set(f, 'Position', [50 60 1050 700], 'Color', 'white');

a1 = axes('Parent', f);
a2 = axes('Parent', f);
a3 = axes('Parent', f);
% a4 = axes('Parent', f);

set(a1, 'Unit', 'normalized', 'Position', [0.06 0.25 0.8 0.2]);
set(a2, 'Unit', 'normalized', 'Position', [0.06 0.49 0.8 0.2]);
set(a3, 'Unit', 'normalized', 'Position', [0.06 0.73 0.8 0.2]);
% set(a4, 'Unit', 'normalized', 'Position', [0.06 0.02 0.8 0.2]);

%initial---
 z = zeros(size_1, size_3);

im1 = image(z,'Parent',a1, 'CDataMapping', 'scaled');
set(a1, 'Clim', clim1);
im2 = image(z,'Parent',a2, 'CDataMapping', 'scaled');
set(a2, 'Clim', clim2);
im3 = image(z,'Parent',a3, 'CDataMapping', 'scaled');
set(a3, 'Clim', clim3);
%---------------
colormap('gray');
N = file_delta;
% size_1 = 128;
% size_3 = 2048;


% fidh = fopen(data_file, 'r');
% h_field = fscanf(fidh,'%e',128*2048*100);
% n_step = length(h_field)/size_1/size_3;
% fclose(fidh);

% length(h_field)

%clim = [-2 2];

movie_filename = strcat(video_path,'field_movie','.avi');
movie_object = VideoWriter(movie_filename);
open(movie_object);

for k = 0:100
  disp(['Processing frame ', num2str(k)]);

  tstart = (k)*N;
  tend = ((k+1)*N-1);
  i = 1;

  fidh1 = fopen([data_file1 num2str(k)], 'r');
  fidh2 = fopen([data_file2 num2str(k)], 'r');
  fidh3 = fopen([data_file3 num2str(k)], 'r');
  h_field1 = fscanf(fidh1,'%e',size_1*size_3*file_delta);
  h_field2 = fscanf(fidh2,'%e',size_1*size_3*file_delta);
  h_field3 = fscanf(fidh3,'%e',size_1*size_3*file_delta);
  %% size(h_field)
  fclose(fidh1);
  fclose(fidh2);
  fclose(fidh3);
  length(h_field1);

  for t = tstart:tend
    local_step = mod(t,file_delta);
    h_field_shot1 = h_field1(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));
    h_field_matrix1 = fliplr(reshape(h_field_shot1,size_3,size_1))';

    h_field_shot2 = h_field2(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));
    h_field_matrix2 = fliplr(reshape(h_field_shot2,size_3,size_1))';

    h_field_shot3 = h_field3(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));
    h_field_matrix3 = fliplr(reshape(h_field_shot3,size_3,size_1))';


    set(im1, 'CData', imfilter(h_field_matrix1,filter));
    set(im2, 'CData', imfilter(h_field_matrix2,filter));
    set(im3, 'CData', imfilter(h_field_matrix3,filter));
    %% imshow(imfilter(h_field_matrix,filter),clim)
    %% figure, surf(h_field_matrix)

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
  %% movie_filename = strcat(video_path,'field_movie_',num2str(t/1e-7,'%3.2f'),'.avi');
  %% movie_object = VideoWriter(movie_filename);
  %% open(movie_object);
  writeVideo(movie_object, mov);
  clear mov;
  %% close(movie_object)
  close(f2);
  %% pack
end

close(movie_object)
