%% Finding the number of rotations within distance range or threshold(hh)
% Ran this script after collection all x1, x2, x3, and x4 values from
% "distances_first.mat" script

%********** anticlockwise rotations **************************************
hh_on = zeros(50,4);hh_off = zeros(50,5);
%finds the first value (between 0 and 0.1)
r1 = nnz(x1>0 & x1<(1/10));rr1 = nnz(y1>0 & y1<(1/10));
r2 = nnz(x2>0 & x2<(1/10));rr2 = nnz(y2>0 & y2<(1/10));
r3 = nnz(x3>0 & x3<(1/10));rr3 = nnz(y3>0 & y3<(1/10));
r4 = nnz(x4>0 & x4<(1/10));rr4 = nnz(y4>0 & y4<(1/10));rr5 = nnz(y5>0 & y5<(1/10));
hh_on(1,:) = [r1 r2 r3 r4];hh_off(1,:) = [rr1 rr2 rr3 rr4 rr5];

%finding the rest of the values for "hh". Range 0.1-0.2/0.2-0.3/...
for u=1:100
    r1 = nnz(x1>(u/10) & x1<(u+1)/10);rr1 = nnz(y1>(u/10) & y1<(u+1)/10);
    r2 = nnz(x2>(u/10) & x2<(u+1)/10);rr2 = nnz(y2>(u/10) & y2<(u+1)/10);
    r3 = nnz(x3>(u/10) & x3<(u+1)/10);rr3 = nnz(y3>(u/10) & y3<(u+1)/10);
    r4 = nnz(x4>(u/10) & x4<(u+1)/10);rr4 = nnz(y4>(u/10) & y4<(u+1)/10);rr5 = nnz(y5>(u/10) & y5<(u+1)/10);
    hh_on(u+1,:) = [r1 r2 r3 r4]; hh_off(u+1,:) = [rr1 rr2 rr3 rr4 rr5];
end

hhh_on = zeros(3,4);hhh_off = zeros(3,5);
b = [1 5;6 10;11 15];           % select first 5 values (0.1-0.5), then the next 5, etc...
for g=1:3
    hhh_on(g,:)= sum(hh_on((b(g,1):b(g,2)),:),1);
    hhh_off(g,:) = sum(hh_off((b(g,1):b(g,2)),:),1);
end

figure   
subplot(3,1,1);bar(hhh_on(1,:))
subplot(3,1,2);bar(hhh_on(2,:))
subplot(3,1,3);bar(hhh_on(3,:))

figure   
subplot(3,1,1);bar(hhh_off(1,:))
subplot(3,1,2);bar(hhh_off(2,:))
subplot(3,1,3);bar(hhh_off(3,:))

%% *********** clockwise rotations ******************************

hh_on = zeros(50,4);hh_off = zeros(50,5);
%finds the first value (between 0 and 0.1)
r1 = nnz(z1>0 & z1<(1/10));rr1 = nnz(zz1>0 & zz1<(1/10));
r2 = nnz(z2>0 & z2<(1/10));rr2 = nnz(zz2>0 & zz2<(1/10));
r3 = nnz(z3>0 & z3<(1/10));rr3 = nnz(zz3>0 & zz3<(1/10));
r4 = nnz(z4>0 & z4<(1/10));rr4 = nnz(zz4>0 & zz4<(1/10));rr5 = nnz(zz5>0 & zz5<(1/10));
hh_on(1,:) = [r1 r2 r3 r4];hh_off(1,:) = [rr1 rr2 rr3 rr4 rr5];

%finding the rest of the values for "hh". Range 0.1-0.2/0.2-0.3/...
for u=1:100
    r1 = nnz(z1>(u/10) & z1<(u+1)/10);rr1 = nnz(zz1>(u/10) & zz1<(u+1)/10);
    r2 = nnz(z2>(u/10) & z2<(u+1)/10);rr2 = nnz(zz2>(u/10) & zz2<(u+1)/10);
    r3 = nnz(z3>(u/10) & z3<(u+1)/10);rr3 = nnz(zz3>(u/10) & zz3<(u+1)/10);
    r4 = nnz(z4>(u/10) & z4<(u+1)/10);rr4 = nnz(zz4>(u/10) & zz4<(u+1)/10);rr5 = nnz(zz5>(u/10) & zz5<(u+1)/10);
    hh_on(u+1,:) = [r1 r2 r3 r4]; hh_off(u+1,:) = [rr1 rr2 rr3 rr4 rr5];
end

hhh_on = zeros(3,4);hhh_off = zeros(3,5);
b = [1 5;6 10;11 15];           % select first 5 values (0.1-0.5), then the next 5, etc...
for g=1:3
    hhh_on(g,:)= sum(hh_on((b(g,1):b(g,2)),:),1);
    hhh_off(g,:) = sum(hh_off((b(g,1):b(g,2)),:),1);
end

figure   
subplot(3,1,1);bar(hhh_on(1,:))
subplot(3,1,2);bar(hhh_on(2,:))
subplot(3,1,3);bar(hhh_on(3,:))

figure   
subplot(3,1,1);bar(hhh_off(1,:))
subplot(3,1,2);bar(hhh_off(2,:))
subplot(3,1,3);bar(hhh_off(3,:))