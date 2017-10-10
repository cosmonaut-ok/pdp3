function inp = m_dct_2d(inp, reverse, dim)

%7/10/2006
%not actually 2d discrete sine transform but a number of 
%1d sine transforms which are applied to all rows/columns of
%input matrix inp (the action of function m_dst_2d is similar
%to the action of fft now)
%reverse = 0, normal dst
%reverse = 1, inverse dst
%dim = 1, dst transform of columns (along y direction)
%dim = 2, dst transform of rows (along x direction)s

[height width] = size(inp);

if dim == 2
    inp = inp';
    [height width] = size(inp);
end

inp = cat(1, inp, flipud(inp(2:end-1,:)));
inp = 1/2*fft(inp,[],1);
inp = real(inp(1:height,:));
if reverse
    inp = inp*2/(height-1);
end

if dim == 2
    inp = inp';
end