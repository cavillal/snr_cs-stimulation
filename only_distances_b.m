%% Calculate distances between each clockwise rotations

%First, load the .mat file containing:
% - clock_rot = index of clockwise rotations/1sec
% - anticlock_rot = inder of anticlock rotations/1sec
% - distance = distance calculated every 1sec (in meters)
% - key = distribution of tests between pre and post injection (1 & 2
% respectively)

% load 'Mus1.mat' 
% Select the desired file to load
filename =uigetfile;
load(filename)

szn = size(distance);               %number of elements in the distance matrix

d = zeros(szn(1),100);
num_rot = zeros(szn(1),1);

%Calculate the distance value for all the turns (d=distances)
for i=1:szn(1)
   a=clock_rot(i,:);                    %select row of clockwise rotation count
   turn = find(a==1);                 %index of where each rotations occurs
   turn1 = [0 turn];                    %adding "0" on the first column for calculating distances
   b = distance(i,:);                   %select each row distance
   c = size(turn);num_rot(i) = c(2);    %number of total rotations per test

    %performs the sum of distances between turns (d = matrix with 
    %all the distances for all the tests)
    for j=1:length(turn1)-1
        s=sum(b((turn1(j)+1:turn1(j+1))));
        d(i,j) = s;
    end
end

d(d==0) = NaN;

separate distance values from pre & post injection using key values matrix 
k = size(key);                          %calculate the size of "key"
pre = zeros(k(1),100);
post = zeros(k(1),100);
    for t=1:length(key)
        pre(t,:) = d(key(t,1),:);
        post(t,:) = d(key(t,2),:);
    end

%eliminate zeros
pre(pre==0) = NaN;
post(post==0) = NaN;

%plot both anti/clock distances pre and post injection
figure
h1 = histogram(pre,'BinWidth',.1);title('anticlockwise rotations')
m = max(max(max(pre),max(post)));
xlabel('distance/rotation');xlim([-.25 m+.25]);
hold
h2 = histogram(post,'BinWidth',.1);
mm1 = max(max(h1.Values),max(h2.Values));
legend('pre','post');
hold off

%Finding the number of rotations within distance range or threshold(hh)
hh = zeros(50,2);
%finds the first value (between 0 and 0.1)
yy = pre>0 & pre<(1/10);yy(yy==0)=[];pre_t=sum(yy);
xx = post>0 & post<(1/10);xx(xx==0)=[];post_t=sum(xx);
hh(1,:) = [pre_t post_t];

%finding the rest of the values for "hh". Range 0.1-0.2/0.2-0.3/...
for u=1:100
    yy = pre>(u/10) & pre<(u+1)/10;yy(yy==0)=[];pre_t=sum(yy);      %finds values between 0.1-0.2 and so on for pre
    xx = post>(u/10) & post<(u+1)/10;xx(xx==0)=[];post_t=sum(xx);   %finds values between 0.1-0.2 and so on for post
    hh(u+1,:) = [pre_t post_t]; 
end

%Plot the number of rotations/distance range
figure
plot(hh);title('clockwise rotations')
legend('pre','post')
xlabel('distanceE-10/rotation')

%% Calculate distances of each anticlockwise rotations
%Calculate the total distances for each turn (d=distances)
d = zeros(szn(1),100);
num_rot = zeros(szn(1),1);
for i=1:szn(1)
   a=anticlock_rot(i,:);                   %select row rotation
   [turn] = find(a==1);                 %where each rotations occurs
   turn1 = [0 turn];                    %adding "0" on the first column
   b = distance(i,:);                   %select row distance
   c = size(turn);num_rot(i) = c(2);    %number of rotations per test

    %performs the sum between turns (d)
    for j=1:length(turn1)-1
        s=sum(b((turn1(j)+1:turn1(j+1))));
        d(i,j) = s;
    end
end

%separate pre & post from usng the values from the key matrix
k = size(key);
pre = zeros(k(1),100);
post = zeros(k(1),100);
for t=1:length(key)
    pre(t,:) = d(key(t,1),:);
    post(t,:) = d(key(t,2),:);
end
%eliminate zeros
pre(pre==0) = NaN;
post(post==0) = NaN;

%plot both anti/clock distances pre and post injection
figure
h1 = histogram(pre,'BinWidth',.1);title('anticlockwise rotations')
m = max(max(max(pre),max(post)));
xlabel('distance/rotation');xlim([-.25 m+.25]);
hold
h2 = histogram(post,'BinWidth',.1);
mm2 = max(max(h1.Values),max(h2.Values));
mm = max(mm1,mm2);
legend('pre','post');ylim([0 mm+2]);
hold off

hh = zeros(50,2);
%finding the first value (between 0 and 0.1)
yy = pre>0 & pre<(1/10);yy(yy==0)=[];pre_t=sum(yy);
xx = post>0 & post<(1/10);xx(xx==0)=[];post_t=sum(xx);
hh(1,:) = [pre_t post_t];

%finding the rest of the values for "hh"
for u=1:100
    yy = pre>(u/10) & pre<(u+1)/10;yy(yy==0)=[];pre_t=sum(yy);      %finds values between 0.1-0.2 and so on
    xx = post>(u/10) & post<(u+1)/10;xx(xx==0)=[];post_t=sum(xx);
    hh(u+1,:) = [pre_t post_t]; 
end
figure
plot(hh);title('anticlockwise rotations')
legend('pre','post')


