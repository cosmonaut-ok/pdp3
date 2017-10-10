% ----------the program shows the plasma electrons flows-------------------
%----for each moment of time the sum of plasma electrons' velocities is----
%------calculated. if sum ~= 0 - the flows are present in the--------------
%---------------------plasma-beam system-----------------------------------
disp ('start');
disp ('Im alive!!');
%--------------files opening for read and write----------------------------
src_id=fopen('c:\dasha\flows\veter-lab-withbeam_s2.txt','r');
dst_id=fopen('D:\Dasha\caverna_new\Levitskij\flows\withbeam_s2_rezults.txt','w');
%---------loop-------------------------------------------------------------
n=0;
tot=0;
while feof(src_id)<1,
    str=fgetl(src_id);
    [q,N]=sscanf(str,'%e');
    M=(N-1)/2;
%     time(n)=q(1);
    Sum=0;
        for i=3:2:N,
            Sum=Sum+q(i);
        end;
%     L=size(V,1);
%         for k=1:L,
%             Sum=Sum+V(k);
%         end;
    Sum=Sum/M;
    tot=tot+Sum;
    disp (q(1));
    disp (Sum);
    fprintf(dst_id,'%e ', q(1));
    fprintf(dst_id,'%e ', Sum);
    fprintf(dst_id,'\n');
    n=n+1;
end;
avtot=tot/n;
disp ('Vt= ');
disp (avtot);
fprintf(dst_id,'DUPA \n');
fprintf(dst_id,'%e ', avtot);
fclose(dst_id); 
fclose(src_id);        
disp ('done!!!');
disp ('DUPA');