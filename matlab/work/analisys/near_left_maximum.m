function var = near_left_maximum(vector, position);

val = vector(position);

i = 1;

while ((position>1)&(vector(position -1) > vector(position)))
    position = position -1;
end

var = position;