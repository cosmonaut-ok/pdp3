function [P p_struct] = beam_injection_1(P,p_struct,N_new_p,dt,vx,ySize,lambda_beam)

x0 = 0;
y0 = ySize*0.48;
lambda = 1e7;

if vx*dt > ySize
    disp('Warning: the beam is too intensive');
    return;
end

free_cells = find(P(3).type == 0);


for k = 1:N_new_p
    rand1 = rand(1);
    rand2 = rand(1);
    P(3).type(free_cells(k)) = 3;
    P(3).x(free_cells(k)) = x0 + vx*dt*rand1;
    P(3).y(free_cells(k)) = y0 + rand2*ySize*0.04;
    P(3).vx(free_cells(k)) = vx;
    P(3).vy(free_cells(k)) = 0;
       
end

p_struct(3).N = length(find(P(3).type));