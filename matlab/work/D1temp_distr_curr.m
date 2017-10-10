%--------------Temperature distribution T(x,t) for 1D model----------------
disp ('start');
disp ('Im alive!!');
%------------------file opening-----------
src_id=fopen('D:\Dasha\caverna_new\Levitskij\Temperatura\odnor_bez_mod_2_ni=ne_s1.txt','r');
dst_id=fopen('D:\Dasha\caverna_new\Levitskij\Temperatura\rezultsS11.txt','w');
% %------------initialize-------------------
%%V_gridt=400;
V_gridx=400; %200 -> 400
xmin=0; xmax=6;
dx=(xmax-xmin)/V_gridx;
position=xmin:dx:xmax;

C1=0.5*9.11*10^-31;
% %------------initial distributions--------
% ------T is temperature distr, K - the quantity of particles in cell------
%for m=1:V_gridx, for n=1:V_gridt, T(m,n)=0; end; end;
%for m=1:V_gridx, for n=1:V_gridt, K(m,n)=0; end; end;%at least one particle in cell
% %-----------------------------------------
% %-----------!!!the main part!!!!!---------
% %-----------------------------------------
n=0;
k=1;
% %--------------LOOP k---------------------------
while feof(src_id)<1,
    str=fgetl(src_id);
    n=n+1;
    if(mod(n,20) ~= 1),
        continue;
    end;
    disp(n);
    [q,N]=sscanf(str,'%e');
    time(k)=q(1);
    %for i=3:2:N,
    %    V((i-1)/2)=q(i); end;
    %for i=2:2:(N-1),
    %    X(i/2)=q(i); end;
    p = 1;
    for i=2:2:N-1,
        X(p) = q(i);
        V(p) = q(i+1);
        p = p + 1;
        %if(mod(p,1000) == 0), disp(p); end;
    end;
    disp('dupa');
    for m=1:(V_gridx-1),
       T(m,k) = 0;
    end;
% %--------------LOOP m---------------------------
    for m=1:(V_gridx-1),
        Vav = 0;
        K = 1;
        for j=1:((N-1)/2),
            if((X(j)>=position(m))&&(X(j)<position(m+1))),
                 K=K+1;
                 Vav=Vav+V(j);
            end;
        end;
        Vav=Vav/K;
%%-------
        for j=1:((N-1)/2),
            if((X(j)>=position(m))&&(X(j)<position(m+1))),
                 T(m,k)=T(m,k)+(V(j)-Vav)^2;
            end;
        end;
        T(m,k)=C1*T(m,k)/(K*10^-18);
    end;
% %--------------END OF THE LOOP m---------------------------
    k=k+1;
end;
% %--------------END OF THE LOOP k---------------------------
V_gridt=k-1;

% %-----------------------------SURF-----------------------------------------
surf(T);
% %------------------------------data writing--------------------------------
fprintf(dst_id,'0 ');
for m=1:(V_gridx-1), fprintf(dst_id,'%e ',position(m));end;
fprintf(dst_id,'\n');
       
for k=1:V_gridt, 
    fprintf(dst_id,'%e ',time(k));
    for m=1:(V_gridx-1), fprintf(dst_id,'%e ',T(m,k));end; 
    fprintf(dst_id,'\n');
end;

% %------------------------cleaning of the variables-------------------------
 %clear T time position X V N K q str;
 disp ('Done!');
fclose(dst_id); 
fclose(src_id);               
disp ('done')