function res = pdp3_movie(raw_data)
mov = [];
% K = ones(399,399,1,12);
for i=1:100
    imshow(squeeze(raw_data(i,:,:)),[-100 100]);
    
%     F(i) = getframe(gca);
            figHandle = gcf;
        frame = getframe(figHandle);
        D(:,:,:,i)   = frame.cdata;
        D(:,:,:,i+1) = frame.cdata;
        D(:,:,:,i+2) = frame.cdata;
        i = i + 3;

%     K(:,:,1,i) = F.cdata(:,:,1);
%   //  mov = addframe(mov,F);
end

f2 = figure('Position', [20 20 100 100]);
    at = axes('Parent', f2);
    mov = immovie(D);
    movie2avi(mov, strcat('E:\taras\_results\','movie4','.avi'));
    clear mov
    close(f2)
    pack

res = 1;