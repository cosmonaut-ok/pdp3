function ret = e_field_intensity(data_file, t_step_num, size_1, size_3, clim);




%clim = [-2 2];
    
        fidh = fopen(data_file, 'r');
        e_field = fscanf(fidh,'%e',size_1*size_3*t_step_num);
        size(e_field);
        fclose(fidh);
        e_field_sum(size_1*size_3)=0;
        e_field_sum =e_field_sum';
    for t = 1:t_step_num
        t
        local_step = mod(t,t_step_num);
        e_field_shot =e_field(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));
        e_field_sum = e_field_sum + e_field_shot.^2;
       
     end
   
 e_field_matrix = fliplr(reshape(e_field_sum,size_3,size_1))';
            
        imshow(e_field_matrix,clim)
      
 ret =  e_field_matrix;
end

