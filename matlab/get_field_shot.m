function res = get_field_shot(path, size_1, size_3, step_t)

fidh = fopen(path, 'r');
h_field = fscanf(fidh,'%e');
n_step = length(h_field)/size_1/size_3;

length(h_field)

h_field_shot = h_field(size_1*size_3*step_t+1:size_1*size_3*(step_t+1));

h_field_matrix = fliplr(reshape(h_field_shot,size_3,size_1))';

res = h_field_matrix;

fclose(fidh);
