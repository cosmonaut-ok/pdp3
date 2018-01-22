function res = get_filed_1d(data_file1,  file_delta, size_1, size_3, r_i);


filter = ones(3,12)/3/12;

path4wr = 'D:\Plasma\_movie\';
name4wr = 'ez';
e_z= zeros(1,size_3);
%N = file_delta;
i=1;
for k = 0:11
 %   file_delta =100;
%    if (k==7)
%        file_delta=95;
 
%    end
%      if (k==9)
%        file_delta=90;
%     end
    
        fidh1 = fopen([data_file1 num2str(k)], 'r');
      
        e_field1 = fscanf(fidh1,'%e',size_1*size_3*file_delta);
     
%         size(h_field)
        fclose(fidh1);
     
        length(e_field1)
    for t = 1:3:(file_delta-1)
        
        %         t        
        local_step = mod(t,file_delta);
        e_field_shot = e_field1(size_1*size_3*local_step+1:size_1*size_3*(local_step+1));

        e_field_matrix = fliplr(reshape(e_field_shot,size_3,size_1))';
        e_field_matrix = medfilt2( e_field_matrix,[23 23]);
        

        e_z (i,1:end)= e_field_matrix(size_1-r_i,1:end);
        i=i+1;
  
        
    end
  


end
 res = e_z;
%res = e_field_matrix;