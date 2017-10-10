function var = near_right_maximum(vector, position);

val = vector(position);

i = 1;

while ((position<length(vector))&(vector(position +1) > vector(position)))
    position = position +1;
end

var = position;