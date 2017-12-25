% function var = rho_movie_create_light3(data_path, video_path, file_delta, size_1, size_3, clim1, clim2, clim3);
function var = rho_movie_create_light3_xml(xml_config_file, video_path, file_delta, clim1, clim2, clim3)
  if ~exist('file_delta', 'var')
    file_delta = 100;
  end

  if ~exist('clim1', 'var')
    clim1 = [0 1];
  end

  if ~exist('clim2', 'var')
    clim2 = [0 1];
  end

  if ~exist('clim3', 'var')
    clim3 = [-1e-7 0];
  end
  
  config_path = fileparts(xml_config_file);
  
  if ~exist('video_path', 'var')
      video_path = config_path;
  end
  
  dom_node = xmlread(xml_config_file);
  dom_root = dom_node.getDocumentElement;
  %% get geometry
  geometry = dom_root.getElementsByTagName('geometry');
  size_1 = str2double(geometry.item(0).getElementsByTagName('n_grid_r').item(0).getFirstChild.getData)-1;
  size_3 = str2double(geometry.item(0).getElementsByTagName('n_grid_z').item(0).getFirstChild.getData)-1;
  
  %% get file to save parameters
  file_save_parameters = dom_root.getElementsByTagName('file_save_paramters');
  local_data_path = file_save_parameters.item(0).getElementsByTagName('path_to_result').item(0).getFirstChild.getData;

  %% disp(class(char(local_data_path)));
  if startsWith(local_data_path, '/')
    data_path = char(local_data_path);
  else
    data_path = [ config_path '/' char(local_data_path) ];
  end
  
  data_file_e1 = strcat(data_path, '/', 'e1');
  data_file_e3 = strcat(data_path, '/', 'e3');
  data_file_rho_beam = strcat(data_path, '/', 'rho_beam');

  filter = ones(3,12); % /3/12; % matrix 12x3, all elements equal to 0.027778

  % create figure window (and place it as current figure)
  f = figure;

  % position: left, bottom, width, height; white background
  set(f, 'Position', [50 60 1050 700], 'Color', 'white');

  %% create axes on figure window f
  a1 = axes('Parent', f);
  a2 = axes('Parent', f);
  a3 = axes('Parent', f);
  
  %% place axes on figure (normalized to whole figure window)
  set(a1, 'Unit', 'normalized', 'Position', [0.06 0.25 0.8 0.2]);
  set(a2, 'Unit', 'normalized', 'Position', [0.06 0.49 0.8 0.2]);
  set(a3, 'Unit', 'normalized', 'Position', [0.06 0.73 0.8 0.2]);

  %% initial---
  z = zeros(size_1, size_3);

  im1 = image(z,'Parent',a1, 'CDataMapping', 'scaled');
  set(a1, 'Clim', clim1);
  im2 = image(z,'Parent',a2, 'CDataMapping', 'scaled');
  set(a2, 'Clim', clim2);
  im3 = image(z,'Parent',a3, 'CDataMapping', 'scaled');
  set(a3, 'Clim', clim3);
  %%---------------
  colormap('gray');

    movie_filename = fullfile(video_path,'field_movie.avi');
    movie_object = VideoWriter(movie_filename);
    open(movie_object);

  N = file_delta; % 100 by default

  for k = 0:100 % TODO: why to 100?
    disp(['Loading files set ', num2str(k)]);

    tstart = k*N; % k=2 -> 200
    tend = ((k+1)*N-1); % k=2 -> 3*100-1 -> 299

      i = 1;
    
    %% Open data files
    fidh_e1 = fopen([data_file_e1 num2str(k)], 'r');
    fidh_e3 = fopen([data_file_e3 num2str(k)], 'r');
    fidh_rho_beam = fopen([data_file_rho_beam num2str(k)], 'r');

    % %e is an exponential notation: ex. 5e2
    h_field_e1 = fscanf(fidh_e1, '%e', size_1*size_3*file_delta); % TODO: why do we times all values?
    h_field_e3 = fscanf(fidh_e3, '%e', size_1*size_3*file_delta);
    h_field_rho_beam = fscanf(fidh_rho_beam, '%e', size_1*size_3*file_delta);

    %% Close data files
    fclose(fidh_e1);
    fclose(fidh_e3);
    fclose(fidh_rho_beam);
    % length(h_field_e1);

    for t = tstart:tend
      local_step = mod(t, file_delta);

      disp(['Processing frame ', num2str(local_step)]);

      h_field_shot1 = h_field_e1(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));
      h_field_matrix1 = fliplr(reshape(h_field_shot1,size_3,size_1))';

      h_field_shot2 = h_field_e3(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));
      h_field_matrix2 = fliplr(reshape(h_field_shot2,size_3,size_1))';

      h_field_shot3 = h_field_rho_beam(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));
      h_field_matrix3 = fliplr(reshape(h_field_shot3,size_3,size_1))';

      set(im1, 'CData', imfilter(h_field_matrix1,filter));
      set(im2, 'CData', imfilter(h_field_matrix2,filter));
      set(im3, 'CData', imfilter(h_field_matrix3,filter));

      %% add colorbar to axes
      colorbar(a1);
      colorbar(a2);
      colorbar(a3);

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
