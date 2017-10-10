function res = temp(inp)

fid = fopen(inp);

time_series = ones(100,399,399);

for i=1:200
  vector = fscanf(fid,'%g',399*399);

  vector = reshape(vector',399,399);

    if i >100

   time_series(i,:,:) = vector';
    end
end

fclose(fid);

res = time_series;