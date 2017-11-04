function var = rho_movie_create_light3(data_path, video_path, file_delta, size_1, size_3, clim1, clim2, clim3);

  if ~exist('file_delta', 'var')
    file_delta = 100;
  end

  if ~exist('size_1', 'var')
    size_1 = 254;
  end

  if ~exist('size_3', 'var')
    size_3 = 2046;
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

  isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

  if isOctave
    try
      pkg load image;
    catch
      pkg install -forge image;
      pkg load image;
    end
  end

  data_file_e1 = strcat(data_path, '/', 'e1');
  data_file_e3 = strcat(data_path, '/', 'e3');
  data_file_rho_beam = strcat(data_path, '/', 'rho_beam');

  filter = ones(3,12)/3/12;

  if isOctave
    f = figure('visible','off');
  else
    f = figure;
  end

  set(f, 'Position', [50 60 1050 700], 'Color', 'white');

  a1 = axes('Parent', f);
  a2 = axes('Parent', f);
  a3 = axes('Parent', f);

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

	if ~isOctave
		movie_filename = strcat(video_path,'field_movie.avi');
		movie_object = VideoWriter(movie_filename);
		open(movie_object);
    end

  N = file_delta;

  for k = 0:100
    disp(['Loading files set ', num2str(k)]);

    tstart = (k)*N;
    tend = ((k+1)*N-1);

		if ~isOctave
			i = 1;
        end

    fidh_e1 = fopen([data_file_e1 num2str(k)], 'r');
    fidh_e3 = fopen([data_file_e3 num2str(k)], 'r');
    fidh_rho_beam = fopen([data_file_rho_beam num2str(k)], 'r');
    %%
    %% %e is an exponential notation
    h_field_e1 = fscanf(fidh_e1, '%e', size_1*size_3*file_delta); % TODO: why do we times all values?
    h_field_e3 = fscanf(fidh_e3, '%e', size_1*size_3*file_delta);
    h_field_rho_beam = fscanf(fidh_rho_beam, '%e', size_1*size_3*file_delta);

    %%
    fclose(fidh_e1);
    fclose(fidh_e3);
    fclose(fidh_rho_beam);
    length(h_field_e1);

    for t = tstart:tend
      local_step = mod(t, file_delta);

      disp(['Processing frame ', num2str(t)]);

      h_field_shot1 = h_field_e1(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));
      h_field_matrix1 = fliplr(reshape(h_field_shot1,size_3,size_1))';

      h_field_shot2 = h_field_e3(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));
      h_field_matrix2 = fliplr(reshape(h_field_shot2,size_3,size_1))';

      h_field_shot3 = h_field_rho_beam(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));
      h_field_matrix3 = fliplr(reshape(h_field_shot3,size_3,size_1))';

      set(im1, 'CData', imfilter(h_field_matrix1,filter));
      set(im2, 'CData', imfilter(h_field_matrix2,filter));
      set(im3, 'CData', imfilter(h_field_matrix3,filter));

      figHandle = gcf;

      %% octave can not write video animation properly
      %% we should just get set of images and then
      %% process it with some other software as video frames
      if isOctave
        img_filename = strcat(video_path, 'field_movie_', int2str(k), '_', int2str(t), '.png');
				print(figHandle, img_filename, '-dpng');
      else
        frame = getframe(figHandle);
				D(:,:,:,i)   = frame.cdata;
				D(:,:,:,i+1) = frame.cdata;
				D(:,:,:,i+2) = frame.cdata;
				i = i + 3;
				drawnow;
      end
    end

		if ~isOctave
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

		if isOctave
			drawnow;
        end

  end
