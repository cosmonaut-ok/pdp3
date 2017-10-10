function res = m_dst(inp, reverse)

[h w] = size(inp);
if (h > 1) 
    inp = inp';
    w = h;
end

inp = cat(2,0,inp,0,-fliplr(inp));
inp = -1/(2*i)*(fft(inp));
res = real(inp(2:w+1));
if (reverse) 
    res = res*2/(w+1);
end

if (h > 1)
    res = res';
end