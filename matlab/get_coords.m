function res = get_coords(path, p_number, step_t)

fidh = fopen(path, 'r');
coords = fscanf(fidh,'%e');
%n_step = length(coords)/size_1/size_3;

length(coords)

coords_shot = coords(p_number*2*step_t+1:p_number*2*(step_t+1));

%h_field_matrix = fliplr(reshape(h_field_shot,size_3,size_1))';

res = coords_shot;

fclose(fidh);
