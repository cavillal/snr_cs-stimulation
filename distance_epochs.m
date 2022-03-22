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

%Calculate the total distances for each turn (d=distances)
for i=1:szn(1)
   a1=anticlock_rot(i,:);a2=clock_rot(i,:);           	%select row rotation
   turn1 = find(a1==1);                                 %index of where anti rotations occurs in a1
   turn2 = find(a2==1);                                 %index of where clock rotations occurs in a2
   b = distance(i,:);                                   %select each row distance
   %number of total rotations ([anticlock clock])
   c = length(turn1);num_rot(i,1) = c;                  %anti_rot
   c = length(turn2);num_rot(i,2) = c;                  %clock_rot
      
   %*********find which turns happen during light ON/light OFF*********
   m = LED(i,:);                    %selec row of light pulses
   m_a = find(m==0);                %find all light off
   m_b = find (m~=0);               %find all light on
 
   mx=m;            %transforming light pulses to binary
   for g=1:length(mx)
       if mx(g)>0
           mx(g)=1;
       end
   end
   
   ff = [];           %finding the cells where light pulses begins and ends
   for g=1:length(mx)-1
       if mx(g+1)-mx(g)~=0 
          ff = [ff g];
       end
   end     
   
   %cells during no light epocs
   nn = [1 ff(1)-1 ff(2)+1 ff(3)-1 ff(4)+1 ff(5)-1 ff(6)+1 ff(7)-1 ff(8)+1 300];
   
   % find anti turns during each light epoch
   l_a1 = turn1((turn1>=ff(1) & turn1<=ff(2)));
   l_a2 = turn1((turn1>=ff(3) & turn1<=ff(4)));
   l_a3 = turn1((turn1>=ff(5) & turn1<=ff(6)));
   l_a4 = turn1((turn1>=ff(7) & turn1<=ff(8)));
   %find clock turns during light epochs
   l_b1 = turn2((turn2>=ff(1) & turn2<=ff(2)));
   l_b2 = turn2((turn2>=ff(3) & turn2<=ff(4)));
   l_b3 = turn2((turn2>=ff(5) & turn2<=ff(6)));
   l_b4 = turn2((turn2>=ff(7) & turn2<=ff(8)));
  
   % find anti turns during no light epochs
   t_a1 = turn1((turn1>=nn(1) & turn1<=nn(2)));
   t_a2 = turn1((turn1>=nn(3) & turn1<=nn(4)));
   t_a3 = turn1((turn1>=nn(5) & turn1<=nn(6)));
   t_a4 = turn1((turn1>=nn(7) & turn1<=nn(8)));
   % find clock turn during no light epoch
   t_b1 = turn2((turn2>=nn(1) & turn2<=nn(2)));
   t_b2 = turn2((turn2>=nn(3) & turn2<=nn(4)));
   t_b3 = turn2((turn2>=nn(5) & turn2<=nn(6)));
   t_b4 = turn2((turn2>=nn(7) & turn2<=nn(8)));
   
   l_anti = {};
   l_clock = {};   
   t_anti = {};
   t_clock = {};
   
   l_anti(i,:) = {l_a1 l_a2 l_a3 l_a4};
   l_clock(i,:) = {l_b1 l_b2 l_b3 l_b4};
   
   t_anti(i,:) = {t_a1 t_a2 t_a3 t_a4};
   t_clock(i,:) = {t_b1 t_b2 t_b3 t_b4};
   
   
end
   

   
   
   
   
   
   
   
   