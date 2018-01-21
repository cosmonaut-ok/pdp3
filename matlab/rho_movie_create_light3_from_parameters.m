function var = rho_movie_create_light3_from_parameters(xml_config_file, clim_e1, clim_e3, clim_rho_beam)

  if ~exist('clim_e1', 'var')
    clim_e1 = [-2e3 2e3];
  end

  if ~exist('clim_e3', 'var')
    clim_e3 = [-2e3 2e3];
  end

  config_path = fileparts(xml_config_file);

  video_path = config_path;
  video_file_name = 'field_movie.avi';

  %% file_delta = 100;

  dom_node = xmlread(xml_config_file);
  dom_root = dom_node.getDocumentElement;
  %% get geometry
  geometry = dom_root.getElementsByTagName('geometry');
  size_1 = str2double(geometry.item(0).getElementsByTagName('n_grid_r').item(0).getFirstChild.getData)-1;
  size_3 = str2double(geometry.item(0).getElementsByTagName('n_grid_z').item(0).getFirstChild.getData)-1;

  y_tick_max = str2double(geometry.item(0).getElementsByTagName('r_size').item(0).getFirstChild.getData);
  x_tick_max = str2double(geometry.item(0).getElementsByTagName('z_size').item(0).getFirstChild.getData);

  %% get normalization parameters
  particles = dom_root.getElementsByTagName('Particles_bunch');
  bunch = dom_root.getElementsByTagName('Particles_bunch');

  %% particles_left_density = str2double(particles.item(0).getElementsByTagName('left_density').item(0).getFirstChild.getData);
  %% particles_right_density = str2double(particles.item(0).getElementsByTagName('right_density').item(0).getFirstChild.getData);
  bunch_density = str2double(particles.item(0).getElementsByTagName('density').item(0).getFirstChild.getData);

  clim_rho_beam = [-(bunch_density*1.6e-19) 0];

  %% get file to save parameters
  file_save_parameters = dom_root.getElementsByTagName('file_save_parameters');
  local_data_path = file_save_parameters.item(0).getElementsByTagName('path_to_result').item(0).getFirstChild.getData;
  frames_per_file = str2num(file_save_parameters.item(0).getElementsByTagName('frames_per_file').item(0).getFirstChild.getData);

  %% disp(class(char(local_data_path)));
  if startsWith(local_data_path, '/')
    data_path = char(local_data_path);
  else
    data_path = [ config_path '/' char(local_data_path) ];
  end

  data_file_e1 = strcat(data_path, '/', 'e1');
  data_file_e3 = strcat(data_path, '/', 'e3');
  data_file_rho_beam = strcat(data_path, '/', 'rho_beam');

  %% filter = ones(3,12)/3/12; % matrix 12x3, all elements equal to 0.027778

  % create figure window (and place it as current figure)
  f = figure('Name', 'PIC Plasma Electron Bunch Modeling Wakefield', 'NumberTitle', 'off');

  % position: left, bottom, width, height; white background
  set(f, 'Position', [50 60 1050 700], 'Color', 'white');

  %% create axes on figure window f
  a1 = axes('Parent', f);
  a2 = axes('Parent', f);
  a3 = axes('Parent', f);

  %% place axes on figure (normalized to whole figure window)
  %% set(a1, 'Unit', 'normalized', 'Position', [0.1 0.74 0.8 0.2]);
  %% set(a2, 'Unit', 'normalized', 'Position', [0.1 0.41 0.8 0.2]);
  %% set(a3, 'Unit', 'normalized', 'Position', [0.1 0.08 0.8 0.2]);
  set(a1, 'Unit', 'normalized', 'Position', [0.1 0.74 0.8 0.2]);
  set(a2, 'Unit', 'normalized', 'Position', [0.1 0.41 0.8 0.2]);
  set(a3, 'Unit', 'normalized', 'Position', [0.1 0.08 0.8 0.2]);

  %%
  %% initial---
  %%

  %%%% initializing images and axes
  z = zeros(size_1, size_3);

  %% set number of marks with names in X and Y axes
  x_ticks_number = 10;
  y_ticks_number = 4;
  %% Titles and Ticks
  x_axe_title = 'Z(m)';
  y_axe_title = 'R(m)';
  colorbar_title = 'V/m';
  x_tick_range = 0:x_tick_max/x_ticks_number:x_tick_max; % we need 10 (or x_ticks_number) ticks
  x_tick_gird_size = 0:size_3/x_ticks_number:size_3;     % from 0 to x_tick_max. it's required
  y_tick_range = 0:y_tick_max/y_ticks_number:y_tick_max; % to convert gird to real size (meters)
  y_tick_gird_size = 0:size_1/y_ticks_number:size_1;     % Same for X and Y axes

  im1 = image(z,'Parent', a1, 'CDataMapping', 'scaled');
  %% set color limits
  a1.CLim = clim_e1;
  a1.CLimMode = 'manual';
  a1.Box = 'off';
  %% set title and axes labels
  a1.Title.String = 'E_r';
  a1.XLabel.String = x_axe_title;
  a1.YLabel.String = y_axe_title;
  a1.TickDir = 'out';
  %% add color bar
  c1 = colorbar(a1);
  c1.TickLabelInterpreter = 'latex';
  c1.Label.String = colorbar_title;
  %% translate gird to real meters in X an Y axes
  a1.XTick = x_tick_gird_size;
  a1.XTickLabel = x_tick_range;
  a1.YTick = y_tick_gird_size;
  a1.YTickLabel = y_tick_range;

  im2 = image(z,'Parent', a2, 'CDataMapping', 'scaled');
  %% set color limits
  a2.CLim = clim_e3;
  a2.CLimMode = 'manual';
  a2.Box = 'off';
  %% set title and axes labels
  a2.Title.String = 'E_z';
  a2.XLabel.String = x_axe_title;
  a2.YLabel.String = y_axe_title;
  a2.TickDir = 'out';
  %% add color bar
  c2 = colorbar(a2);
  c2.TickLabelInterpreter = 'latex';
  c2.Label.String = colorbar_title;
  %% translate gird to real meters in X an Y axes
  a2.XTick = x_tick_gird_size;
  a2.XTickLabel = x_tick_range;
  a2.YTick = y_tick_gird_size;
  a2.YTickLabel = y_tick_range;

  im3 = image(z,'Parent', a3, 'CDataMapping', 'scaled');
  %% set color limits
  a3.CLim = clim_rho_beam;
  a3.CLimMode = 'manual';
  a3.Box = 'off';
  %% set title and axes labels
  a3.Title.String = 'Electrons density';
  a3.XLabel.String = x_axe_title;
  a3.YLabel.String = y_axe_title;
  a3.TickDir = 'out';
  %% add color bar
  c3 = colorbar(a3);
  c3.TickLabelInterpreter = 'latex';
  c3.Label.String = 'm^{-3}';
   % translate e_bunch to electrons density and set custom tick labels
   % and direction (because electron charge is negative). Also, round
   % tick labels array to nearest rounded values for short marks
  c3.TickLabels = round(0:bunch_density/3:bunch_density, 2, 'significant');
  c3.Direction = 'reverse';
  %% translate gird to real meters in X an Y axes
  a3.XTick = x_tick_gird_size;
  a3.XTickLabel = x_tick_range;
  a3.YTick = y_tick_gird_size;
  a3.YTickLabel = y_tick_range;

  %%---------------
  colormap('gray');

  movie_filename = fullfile(video_path,video_file_name);
  movie_object = VideoWriter(movie_filename);
  open(movie_object);

  N = frames_per_file; % 100 by default

  for k = 0:100 % TODO: why to 100?

    tstart = k*N; % k=2 -> 200
    tend = ((k+1)*N-1); % k=2 -> 3*100-1 -> 299

    i = 1;

    if ~(exist([data_file_e1 num2str(k)], 'file') == 2) || ~(exist([data_file_e3 num2str(k)], 'file') == 2) || ~(exist([data_file_rho_beam num2str(k)], 'file') == 2)
      disp('No more data files exists. Exiting');
      return
    end

    disp(['Loading files set ', num2str(k)]);

    %% Open data files
    fidh_e1 = fopen([data_file_e1 num2str(k)], 'r');
    fidh_e3 = fopen([data_file_e3 num2str(k)], 'r');
    fidh_rho_beam = fopen([data_file_rho_beam num2str(k)], 'r');

    % %e is an exponential notation: ex. 5e2
    h_field_e1 = fscanf(fidh_e1, '%e', size_1*size_3*frames_per_file); % TODO: why do we times all values?
    h_field_e3 = fscanf(fidh_e3, '%e', size_1*size_3*frames_per_file);
    h_field_rho_beam = fscanf(fidh_rho_beam, '%e', size_1*size_3*frames_per_file);

    %% Close data files
    fclose(fidh_e1);
    fclose(fidh_e3);
    fclose(fidh_rho_beam);
    % length(h_field_e1);

    for t = tstart:tend
      local_step = mod(t, frames_per_file);

      disp(['Processing frame ', num2str(local_step)]);

      h_field_shot1 = h_field_e1(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));
      h_field_matrix1 = fliplr(reshape(h_field_shot1,size_3,size_1))';

      h_field_shot2 = h_field_e3(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));
      h_field_matrix2 = fliplr(reshape(h_field_shot2,size_3,size_1))';

      h_field_shot3 = h_field_rho_beam(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));
      h_field_matrix3 = fliplr(reshape(h_field_shot3,size_3,size_1))';

      %% set(im1, 'CData', imfilter(h_field_matrix1,filter));
      %% set(im2, 'CData', imfilter(h_field_matrix2,filter));
      %% set(im3, 'CData', imfilter(h_field_matrix3,filter));

      set(im1, 'CData', h_field_matrix1);
      set(im2, 'CData', h_field_matrix2);
      set(im3, 'CData', h_field_matrix3);

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
