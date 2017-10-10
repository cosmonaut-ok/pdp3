src_id=fopen('D:\Dasha\caverna_new\Levitskij\Temperatura\rezultsS1.txt','r');
disp ('begin');
str=fgetl(src_id);
[q,N]=sscanf(str,'%e');
for i=2:N, x(i-1)=q(i); end;
%disp('dimx=');
%disp(length(x));
%disp ('end');
j=1;
while feof(src_id)<1,
     str=fgetl(src_id);
     [q,N]=sscanf(str,'%e'); 
     t(j)=q(1);
     for i=2:N, F(j,(i-1))=q(i); end;
     j=j+1;
end;
disp('dimt=');
disp(length(t));
disp('dimx=');
disp(length(x));
disp('dimF=');
disp(size(F));
fclose(src_id);
surf(x,t,F);
% clear x y dW N X Y VX VY q str;
     