%% Manual mouse head tracker

fr = 5433;              %first frame light pulse R1:1021-1950;R2:2780-3659;R3:4536-5422;R4:6298-7190
fr_end = 5850;          %last frame light pulse
g = fr_end-fr;          %number of frames
hd_x = [];hd_y = [];
light = [];
dd = [];
angle = [];

v = VideoReader('Test_6.mp4');

for i = 0:g
    %frame = read(v,(fr-20)+i); reads it from 20 frames before first light
    frame = read(v,fr+i);       %reads it from the first light pulse
    image(frame)
    [x,y,button] = ginput(1);
    hd_x = [hd_x x];            %x-axis coordinates
    hd_y = [hd_y y];            %y-axis coordinates
    light = [light button];     %1=(left button):no light - 3=(right button):light
end

plot(hd_x,hd_y)                 %plot tracked head movements

%Calculating distance between points
for k=1:length(hd_x)
    d = sqrt((hd_x(k+1)-hd_x(k))^2+(hd_y(k+1)-hd_y(k))^2);
    dd = [dd d];
end

%Calculating angle between vector
for h =1:length(hd_x)
    u = [hd_x(h) hd_x(h+1) hd_x(h+2)];
    v = [hd_y(h) hd_y(h+1) hd_y(h+2)];
    ThetaInDegrees = atan2d(norm(cross(u,v)),dot(u,v));
    angle = [angle ThetaInDegrees];
end

%eliminating the last 2 elements of light (matching matrices)
if length(light)-length(angle) ~= 0  
    j = length(light)-length(angle);
    for jj = 1:j
        light(end)=[];
    end
end

%ploting head turn distance
plot (dd);
yyaxis right
plot(light)

%plotting head turn angle
figure
plot(angle)
yyaxis right
plot(light)

l_off = angle(light == 1);
l_on = angle(light == 3);
d_off = dd(light == 1);
d_on = dd(light ==3);

save ('11-26-19-test4_R1.mat','hd_x','hd_y','light','l_off','l_on','d_off','d_on')

%d = sqrt((x2-x1)^2+(y2-y1)^2); distance formula