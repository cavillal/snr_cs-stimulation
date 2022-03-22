close all
g = find (animals == 7);
d = distance_5;
d(d==0) = 0.000000001;
abs_clock = clock_angle_5(g,:)./d(g,:);
abs_anticlock = anticlock_angle_5(g,:)./d(g,:);
plot (abs_clock)
hold
plot (abs_anticlock)

bin_clock = max(abs_clock)/10;
bin_anticlock = max(abs_anticlock)/10;
thres_clock = zeros(1,10);
thres_anticlock = zeros(1,10);
for i=4
    m = abs_clock > i*1e11;
    m1 = sum(m);
    n = abs_anticlock > i*1e11;
    n1 = sum(n);
    thres_clock = [thres_clock m1]
    thres_anticlock = [thres_clock n1]
end
