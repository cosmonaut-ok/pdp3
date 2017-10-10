function res = get_filed_1d(data_file1,  file_delta, size_1, size_3, r_i);




path4wr = 'D:\Plasma\_movie\';
name4wr = 'ez';
e_z= zeros(1,size_3);
N = file_delta;
for k = 0:1

    i = 1;

        fidh1 = fopen([data_file1 num2str(k)], 'r');
      
        e_field1 = fscanf(fidh1,'%e',size_1*size_3*file_delta);
     
%         size(h_field)
        fclose(fidh1);
     
        length(e_field1)
    for t = 1:(N-1) 
        %         t
        n_st=(r_i-1)*(size_3)+1;
        n_end = (r_i)*(size_3);
        ez_t = e_field1((n_st+size_1*size_3*t):(n_end+size_1*size_3*t));
        ez_t = medfilt2(ez_t,20);
        e_z (t,1:end)=ez_t;
  
        
    end
  


end
 res = e_z;
