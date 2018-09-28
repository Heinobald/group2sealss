function [output1,output2, output3] = smoothday(variable, depth, time)
%   Builds mean of variable for every day 
%   variable(depth, time)
%   depth
%   d = datenum format 
%   output1 = smoothed variable
%   output2 = smoothed time (as array)
temp=transpose(variable); %time needs to be in y direction
d=datetime(time, 'ConvertFrom', 'datenum'); %timetable function needs datetime format

for t=1:length(depth)
   tt=timetable(d,temp(:,t));
    tt=retime(tt,'daily', 'mean');  
    tt=timetable2table(tt); %row 1: time row 2: salinity at depth t
    timen=tt(:,1); %new date format showing each day with data once
    tt(:,1)=[]; %delete time column out of table
    tt=table2array(tt); %transform sal table (1x:) to sal array (1x:)
    output1(t,:)=tt; %add sal column to a new array sal47s2
end  
output3=table2array(timen);
output2=datenum(output3); %gives the time as datenum
end

