%Heiner
%13.09.2018 excercise statistics

[file, txt, raw]=xlsread('C:\Users\Heiner\Desktop\Marine Project\0913_power_analysis\Lindegarth_lectures and practicals\bqi_data.xlsx');

    % years=file(:,1);
    % y2002=file(years==2002,:); %only data for year 2002
y2002=raw(file(:,1)==2002,:); %reduce the 2 rows above into one
y2003=file(file(:,1)==2003,:);
test=cell2mat(raw(2:161,2));
test=str2mat(y2002(:,6),4);