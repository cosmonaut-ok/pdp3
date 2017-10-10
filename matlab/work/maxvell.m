V_grid=50;
src_id=fopen('C:\dasha\odnor_bez_mod_2_ni=ne_s1.txt','r');
% dst_id=fopen('I:\Anisimov\Temperature\distribution_dyn.txt','w');
%--------------------------------------------------------------------------
disp ('I am alive');
str=fgetl(src_id);
[q,N]=sscanf(str,'%e');
V_min=q(3);
V_max=q(3);
i=3; while i<N,if V_min>q(i),V_min=q(i);end; if V_max<q(i),V_max=q(i);end;i=i+2;end;
while feof(src_id)<1,
   i=3;
   str=fgetl(src_id);
  [q,N]=sscanf(str,'%e');
   while i<N,if V_min>q(i),V_min=q(i);end; if V_max<q(i),V_max=q(i);end;i=i+2;end;
end;
dV=(V_max-V_min)/V_grid;
%%--------------------------------------------------------------------------
fseek(src_id,0,-1);
%--------------------------------------------------------------------------
t=1;
while feof(src_id)<1,
   str=fgetl(src_id);
   [q,N]=sscanf(str,'%e');
   for j=1:V_grid,F_v(t,j)=0;end;
   i=3;   
   while i<N,
      for j=1:V_grid,if (q(i)>=V_min+dV*j)&&(q(i)<V_min+dV*(j+1)), F_v(t,j)=(F_v(t,j)+1);end;end;
      i=i+2;
   end;
   t=t+1;
%    V=V_min:dV:V_max-dV;
%    for j=1:V_grid,fprintf(dst_id,'%e %e %e ',q(1),V(j),F_v(j)/(0.5*N));end; 
%    fprintf(dst_id,'\n');
%--------------------------------------------------------------------------
end;
disp (V_max);
disp (V_min);
surf(F_v);
clear Fn  q str;
disp ('Done!');
% fclose(dst_id); 
fclose(src_id);
%--------------------------------------------------------------------------
% src_id=fopen(I:\Anisimov\Temperature\distribution_dyn.txt','r');
% j=1;
% while feof(src_id)<1,
%  str=fgetl(src_id);
%   [q,N]=sscanf(str,'%e');
%   i=1;
%   for l=1:3:N,
%       T(j)=q(l);
%       V(i)=q(l+1);
%       F_vt(j,i)=q(l+2);
%       i=i+1;
%  end;
%   j=j+1;
% end;
% fclose(src_id);
% surf(V,T,F_vt);
% clear T F_vt V q str;
% 
