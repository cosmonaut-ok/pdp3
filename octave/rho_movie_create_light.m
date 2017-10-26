function var = rho_movie_create_light(data_file, file_delta, size_1, size_3, clim);

w = 1.0e9;
Tmod = 2*pi/w;

filter = ones(9,12)/9/12;




path4wr = 'e:\Science[Plasma]\_movie\';
name4wr = 'h2';



t_begin = 1.6e-6;
dt = 1e-9;
var = 'rho_sp';

f = figure;

N = file_delta;
% size_1 = 128;
% size_3 = 2048;


% fidh = fopen(data_file, 'r');
% h_field = fscanf(fidh,'%e',128*2048*100);
% n_step = length(h_field)/size_1/size_3;
% fclose(fidh);

% length(h_field)


%clim = [-2 2];

for k = 1:100

    tstart = (k-1)*N;
    tend = (k*N-1);
    i = 1;

        fidh = fopen([data_file num2str(k-1)], 'r');
        h_field = fscanf(fidh,'%e',size_1*size_3*file_delta);
        size(h_field)
        fclose(fidh);
    for t = tstart:tend
        t
        local_step = mod(t,file_delta);
        h_field_shot = h_field(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));

        h_field_matrix = fliplr(reshape(h_field_shot,size_3,size_1))';
            
        imshow(h_field_matrix,clim)
       % figure, surf(h_field_matrix)



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
    movie2avi(mov, strcat(path4wr,name4wr,'_',num2str(t/1e-7,'%3.2f'),'.avi'));
    clear mov
    close(f2)
    pack



end
