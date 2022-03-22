%% Analysis of the distances/frame during, between and before light pulses
%l = head angle difference between frames
%d = had distance between frames

file = uigetfile;
load(file)

szn = size(hd_x);
l_on(l_on==0) = [];
l_off(l_off==0) = [];
d_on(d_on==0) = [];
d_off(d_off==0) = [];

%Calculate means
m_l_off = mean(l_off);
m_l_on = mean(l_on);
m_d_off = mean(d_off);
m_d_on = mean(d_on);

figure
subplot(4,2,[1 2])
histogram(l_off,'BinWidth',.5);hold
histogram(l_on,'BinWidth',.5); hold off
legend('ligh OFF','light ON');title(file);
subplot(4,2,[3 4])
histogram(d_off,'BinWidth',.5);hold
histogram(d_on,'BinWidth',.5); hold off
legend('ligh OFF','light ON')

%Calculate standard deviations
s_l_off = std(l_off);
s_l_on = std(l_on);
s_d_off = std(d_off);
s_d_on = std(d_on);

%% Bin events from angle 
total_events = [length(l_off) length(l_on)];        %all events together (ON/OFF)
subplot(4,2,5)
bar (total_events)
xlabel(['OFF' 'ON'])

%bin events for angle (from 0 to 1SD)
xx = zeros(1,10);yy = zeros(1,10);
xl= l_off < s_l_off;xl(xl==0) = [];
yl= l_on < s_l_on;yl(yl==0) = [];
xx(1) = length(xl);
yy(1) = length(yl);

%bin events (from 1SD and so on)
for h=1:9
    xl = l_off > s_l_off*h & l_off < s_l_off*h+1;xl(xl==0)=[];
    yl = l_on > s_l_on*h & l_on < s_l_on*h+1;yl(yl==0)=[];
    
    xx(1,h+1) = length(xl);
    yy(1,h+1) = length(yl);
end 
    
perc = zeros(1,10);
for g=1:10
    perc_off = xx(g)/(xx(g)+yy(g));perc_on = 1-perc_off;
    perc(1,g) = perc_off;
end

subplot (4,2,7)
bar (perc);
ylim([.5 1])

%% Bin events from distance
total_events = [length(d_off) length(d_on)];        %all events together (ON/OFF)
subplot(4,2,6)
bar (total_events)

% bin events for distance (from 0 to 1SD)
xx = zeros(1,10);yy = zeros(1,10);
xd= d_off < s_d_off;xd(xd==0) = [];
yd= d_on < s_d_on;yd(yd==0) = [];
xx(1) = length(xd);
yy(1) = length(yd);

%bin events (from 1SD and so on)
for h=1:9
    xd = d_off > s_d_off*h & d_off < s_d_off*h+1;xd(xd==0)=[];
    yd = d_on > s_d_on*h & d_on < s_d_on*h+1;yd(yd==0)=[];
    
    xx(1,h+1) = length(xd);
    yy(1,h+1) = length(yd);
end 
    
perc = zeros(1,10);
for g=1:10
    perc_off = xx(g)/(xx(g)+yy(g));perc_on = 1-perc_off;
    perc(1,g) = perc_off;
end

subplot (4,2,8)
bar (perc);
ylim([.5 1])





