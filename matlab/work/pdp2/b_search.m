function res = b_search(x,vector)

% returns the index of vector element which
% is the closest tho the number x 

N = length(vector);

switch N
    case 0
        res = 0;
        disp('The vector is empty')
        return
        
    case 1
        res = 1;
        disp('Scalar')
        return
        
    otherwise

        if x > vector(end)
            res = N;
            return
        end
        if x < vector(1)
            res = 1;
            return
        end
        l = 1;
        r = N;
        while r - l > 1
            mid = floor((r+l)/2);
            if x < vector(mid)
                r = mid;
            else
                l = mid;
            end
        end
                      
end
% if abs(x - vector(r)) > abs(x - vector(l))
%     res = l;
% else
%     res = r;
% end
res = l;
