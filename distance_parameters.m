%% This file adds the distances values within and outside light pulses -
%it uses the previously calculated values from the analysis manual_head_tracker 
file = uigetfile;
load(file)

plot(hd_x,hd_y)                 %plot tracked head movements

% dd = [];
% 
% %Calculating distance between points
% for k=1:length(hd_x)-1
%     d = sqrt((hd_x(k+1)-hd_x(k))^2+(hd_y(k+1)-hd_y(k))^2);
%     dd = [dd d];
% end
% 
% %eliminating the last 2 elements of light (matching matrices)
% if length(light)-length(dd) ~= 0  
%     j = length(light)-length(dd);
%     for jj = 1:j
%         light(end)=[];
%     end
% end
% 
% d_off = dd(light == 1);
% d_on = dd(light ==3);
% 
% save (file,'hd_x','hd_y','light','l_off','l_on','d_off','d_on')
% 
% clear all