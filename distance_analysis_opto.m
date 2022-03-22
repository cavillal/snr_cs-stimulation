%% Calculate distances between each anti and clockwise rotations during light on and off

%First, load the .mat file containing:
% - clock_rot = index of clockwise rotations/1sec
% - anticlock_rot = inder of anticlock rotations/1sec
% - distance = distance calculated every 1sec (in meters)
% - LED = matrix where light was on/off every 1sec

filename = uigetfile;
load(filename)
savename = append("data_",filename);

szn = size(distance);               %number of elements in the distance matrix

d1 = zeros(szn(1),100);d2=d1;
num_rot = zeros(szn(1),2);
turn1 = zeros(10,1);turn2=turn1;
total_turns = zeros(szn(1),4);
% session = 7;                         %number of the column for the session to analyze

%Calculate the total distances for each turn (d=distances)
for i=1:szn(1)
   a1=anticlock_rot(i,:);a2=clock_rot(i,:);           %select row rotation
   turn1 = find(a1==1);                 %index of where anti rotations occurs in a1
   turn2 = find(a2==1);                 %index of where clock rotations occurs in a2
   b = distance(i,:);                   %select each row distance
   %number of total rotations ([anticlock clock])
   c = length(turn1);num_rot(i,1) = c;      %anti_rot
   c = length(turn2);num_rot(i,2) = c;      %clock_rot
     
   
   %*********find which turns happen during light ON/light OFF*********
   m = LED(i,:);                    %selec row of light pulses
   m_a = find(m==0);                %find all light off
   m_b = find (m~=0);               %find all light on
 
%    mx=m;            %transforming light pulses to binary
%    for g=1:length(mx)
%        if mx(g)>0
%            mx(g)=1;
%        end
%    end
%    
%    ff = [];           %finding the cells where light pulses begins and ends
%    for g=1:length(mx)-1
%        if mx(g+1)-mx(g)~=0 
%           ff = [ff g];
%        end
%    end     
%    
%    % find anti turns during each light epoch
%    t_b1 = turn1((turn1>=ff(1) & turn1<=ff(2)));
%    t_b2 = turn1((turn1>=ff(3) & turn1<=ff(4)));
%    t_b3 = turn1((turn1>=ff(5) & turn1<=ff(6)));
%    t_b4 = turn1((turn1>=ff(7) & turn1<=ff(8)));
%    
%    tt = {};
%    tt (i,:) = {t_b1 t_b2 t_b3 t_b4};
   
   anti_l_on = zeros(1,20);
   anti_l_off = zeros(1,20);
   clock_l_on = zeros(1,20);
   clock_l_off = zeros(1,20);
   
   for h=1:length(turn1)            % for anticlock rotations
       t_b = find(turn1(h)==m_b);   %find anticlock rot in light on
       t_a = find(turn1(h)==m_a);   %find anticlock rot in light off       
       if isempty(t_b)==1
           t_b=0;
       end
       if isempty(t_a)==1
           t_a = 0;
       end
       anti_l_on(h) = t_b;anti_l_off(h) = t_a; %[on off]
   end
   
   for h=1:length(turn2)            % for clock rotations
       t_b = find(turn2(h)==m_b);   %find anticlock rot in light on
       t_a = find(turn2(h)==m_a);   %find anticlock rot in light off       
       if isempty(t_b)==1
           t_b=0;
       end
       if isempty(t_a)==1
           t_a = 0;
       end
       clock_l_on(h) = t_b;clock_l_off(h) = t_a; %[on off]
   end
   
   anti_l_on(anti_l_on==0)=[];
   anti_l_off(anti_l_off==0)=[];
   clock_l_on(clock_l_on==0)=[];
   clock_l_off(clock_l_off==0)=[];  
 
   turn_anti_off = zeros(1);
   turn_anti_on = zeros(1);
   turn_clock_off = zeros(1);
   turn_clock_on = zeros(1);
 
   
   for h=1:length(anti_l_on)
       xx = m_b(anti_l_on(1,h));
       turn_anti_on(h) = xx;
   end
   for h=1:length(anti_l_off)
       xx = m_a(anti_l_off(1,h));
       turn_anti_off(h) = xx;
   end
   for h=1:length(clock_l_on)
       xx = m_b(clock_l_on(1,h));
       turn_clock_on(h) = xx;
   end
   for h=1:length(clock_l_off)
       xx = m_a(clock_l_off(1,h));
       turn_clock_off(h) = xx;
   end
    
   %quantification of all the turns (anti on;anti off;clock on;clockoff)
   tt = [nnz(turn_anti_on) nnz(turn_anti_off) nnz(turn_clock_on) nnz(turn_clock_off)];
   total_turns(i,:) = tt;
 
   %adding a 0 in the first column for all the turnings
   [turn_anti_off] = [0 turn_anti_off];
   [turn_anti_on] = [0 turn_anti_on];
   [turn_clock_off] = [0 turn_clock_off];
   [turn_clock_on] = [0 turn_clock_on];
    
   d_anti_off = (zeros(1,length(anti_l_off)));
   d_anti_on = (zeros(1,length(anti_l_on)));
   d_clock_off = (zeros(1,length(clock_l_off)));
   d_clock_on = (zeros(1,length(clock_l_on)));
   
   
    %performs the sum of distances between turns (d = matrix  with 
    %all the distances for all the tests)
    for j=1:length(turn_anti_off)-1
        s=sum(b((turn_anti_off(j)+1:turn_anti_off(j+1))));
        d_anti_off(j) = s;
    end
    
    for j=1:length(turn_anti_on)-1
        s=sum(b((turn_anti_on(j)+1:turn_anti_on(j+1))));
        d_anti_on(j) = s;
    end
    
    for j=1:length(turn_clock_off)-1
        s=sum(b((turn_clock_off(j)+1:turn_clock_off(j+1))));
        d_clock_off(j) = s;
    end
    
    for j=1:length(turn_clock_on)-1
        s=sum(b((turn_clock_on(j)+1:turn_clock_on(j+1))));
        d_clock_on(j) = s;
    end
    
    %collection of all the distances
    %data = columns - tests. rows 1:anti_off 2:anti_on 3:clock_off 4:clock_on
    data (i,:) = {d_anti_on d_anti_off d_clock_on d_clock_off}; 
    
end

   xname = ["anti-on" "anti-off" "clock-on" "clock-off"];
   bar(total_turns');set(gca,'xticklabel',xname);
   title(filename)

% save data filename and the total number of turns between anti and clock
% turns
save (savename,'data','filename','total_turns')   %saves the data into savename

clear all

%%


% keep total_turns d_anti_off d_anti_on d_clock_off d_clock_on
% save ('101719','filename','session','total_turns','d_anti_off','d_anti_on','d_clock_off','d_clock_on')

% 
% d1(d1==0) = NaN;
% d2(d2==0) = NaN;

% %separate distance values from pre & post injection using key values matrix 
% k = size(key);                          %calculate the size of "key"
% pre = zeros(k(1),100);
% post = zeros(k(1),100);
% for t=1:length(key)
%     pre(t,:) = d(key(t,1),:);
%     post(t,:) = d(key(t,2),:);
% end
% 
% %eliminate zeros
% pre(pre==0) = NaN;
% post(post==0) = NaN;
% 
% %plot both anti/clock distances pre and post injection
% figure
% h1 = histogram(pre,'BinWidth',.1);title('anticlockwise rotations')
% m = max(max(max(pre),max(post)));
% xlabel('distance/rotation');xlim([-.25 m+.25]);
% hold
% h2 = histogram(post,'BinWidth',.1);
% mm1 = max(max(h1.Values),max(h2.Values));
% legend('pre','post');
% hold off
% 
% %Finding the number of rotations within distance range or threshold(hh)
% hh = zeros(50,2);
% %finds the first value (between 0 and 0.1)
% yy = pre>0 & pre<(1/10);yy(yy==0)=[];pre_t=sum(yy);
% xx = post>0 & post<(1/10);xx(xx==0)=[];post_t=sum(xx);
% hh(1,:) = [pre_t post_t];
% 
% %finding the rest of the values for "hh". Range 0.1-0.2/0.2-0.3/...
% for u=1:100
%     yy = pre>(u/10) & pre<(u+1)/10;yy(yy==0)=[];pre_t=sum(yy);      %finds values between 0.1-0.2 and so on for pre
%     xx = post>(u/10) & post<(u+1)/10;xx(xx==0)=[];post_t=sum(xx);   %finds values between 0.1-0.2 and so on for post
%     hh(u+1,:) = [pre_t post_t]; 
% end
% 
% %Plot the number of rotations/distance range
% figure
% plot(hh);title('clockwise rotations')
% legend('pre','post')
% xlabel('distanceE-10/rotation')

