szn = size(distance);

% %plot all the distances
% for i=1:szn(1)
%     a=distance(i,:);
%     figure
%     plot(a)
% end

%# of rotations
for i=1:szn(1)
    sum(clock_rot(i,:))
    sum(anticlock_rot(i,:))
    %pause
end

% Calculation of distance between turns
d=[];
%d = zeros(szn(1),length(turn));
for i=1:szn(1)
   a=clock_rot(i,:);                   %select row rotation
   [turn] = find(a==1);
   turn1 = [0 turn];
   b = distance(i,:);                   %select row distance

    %performs the sum between turns (d)
    for j=1:length(turn1)-1
        s=sum(b((turn1(j)+1:turn1(j+1))));
        d(i,j) = s;
    end
end

plot(d)

yy=[];
%cummulative 
for h=1:szn(1)
    for g=1:szn(1)
        y=sum(d(h,:)<g/10 & d(h,:)>0.00001);
    yy(h,g) = y;
    end
end

for h=1:szn(1)
    for g=1:szn(1)
        y=sum(d(h,:)<g/20 & d(h,:)>(g/20)-0.05);
    yy(h,g) = y;
    end
end


y = sum(d(1,:)>5/10);                   %find number of d > 0.4


for g=1:25
    y=sum(d(10,:)<g/10 & );
    yy(1,g) = y;
end


    
    
 
