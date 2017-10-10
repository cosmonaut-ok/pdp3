function var = iteration_vs_superscalar;

N = 10000000;


t = rand(1,N);
t1 = zeros(1,N);


tic

t1 = exp(t.*t) + 2;

disp('superscalar time:')
time1 = toc

z = rand(1,N);


tic
for i = 1:N
    t1(i) = exp(z(i)*z(i)) + 3;
end
disp('iteration time:');
time2 = toc

var = time1 - time2;
