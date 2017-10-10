function res = m_dct(inp, reverse)

[h w] = size(inp);
if (h > 1) 
    inp = inp';
    w = h;
end

inp = cat(2,inp,fliplr(inp(2:end-1)));
inp = 1/2*fft(inp);
res = real(inp(1:w));
if (reverse) 
    res = res*2/(w-1);
end

if (h > 1)
    res = res';
end