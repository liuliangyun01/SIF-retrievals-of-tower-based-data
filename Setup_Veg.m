% ------设置初始植被高度，最终植被高度，生长期初始时间-------- -% 
% ---------为了获得上下行路径(H_tower-H_veg)--------- -% 
function [Time_init, Total_days, Veg_0, Veg_max] = Setup_Veg(Name_Tower, Year, Date)
    Vegetations = [];
    if strcmp(Name_Tower, 'XTS')   
        if str2num(Year) == 2017
            if (20170429 <= str2num(Date)) && (str2num(Date) <= 20170608)  % 小麦
                Time_init = datenum(2017, 04, 29); 
                Total_days = datenum(2017, 06, 08) - datenum(2017, 04, 29);       
                Veg_0 = 0.7; 
                Veg_max = 0.85; 
            elseif (20170706 <= str2num(Date)) && (str2num(Date) <= 20170912)  % 玉米
                Time_init = datenum(2017, 07, 06); 
                Total_days = datenum(2017, 09, 12) - datenum(2017, 07, 06);   
                Veg_0 = 0; 
                Veg_max = 2.5; 
            elseif (20171028 <= str2num(Date)) && (str2num(Date) <= 20171231)  % 小麦
                Time_init = datenum(2017, 10, 28); 
                Total_days = datenum(2017, 12, 31) - datenum(2017, 10, 28);   
                Veg_0 = 0; 
                Veg_max = 0.1; 
            else
                Time_init = 0;
                Total_days = 0;
                Veg_0 = 0; 
                Veg_max = 0;
            end
        elseif str2num(Year) == 2018
            if (20180101 <= str2num(Date)) && (str2num(Date) <= 20180615)  % 小麦
                Time_init = datenum(2018, 01, 01); 
                Total_days = datenum(2018, 06, 15) - datenum(2018, 01, 01)+ 1;      % +1是为了和下面的0，区分开，为了将来判断用
                Veg_0 = 0.1; 
                Veg_max = 0.85; 
            elseif (20180628 <= str2num(Date)) && (str2num(Date) <= 20180921)  % 玉米
                Time_init = datenum(2018, 06, 28); 
                Total_days = datenum(2018, 09, 21) - datenum(2018, 06, 28)+ 1;  
                Veg_0 = 0; 
                Veg_max = 2.5; 
            elseif (20181024 <= str2num(Date)) && (str2num(Date) <= 20181231)  % 小麦
                Time_init = datenum(2018, 10, 24); 
                Total_days = datenum(2018, 12, 31) - datenum(2018, 10, 24)+ 1;  
                Veg_0 = 0; 
                Veg_max = 0.1; 
            else
                Time_init = 0;
                Total_days = 0;  
                Veg_0 = 0; 
                Veg_max = 0;
            end
        elseif str2num(Year) == 2019
            if (20190101 <= str2num(Date)) && (str2num(Date) <= 20190611)  % 小麦
                Time_init = datenum(2019, 01, 01); 
                Total_days = datenum(2019, 06, 11) - datenum(2019, 01, 01)+ 1;  
                Veg_0 = 0.1; 
                Veg_max = 0.85; 
            elseif (20190816 <= str2num(Date)) && (str2num(Date) <= 20190915)  % 杂草
                Time_init = datenum(2019, 08, 16); 
                Total_days = datenum(2019, 09, 15) - datenum(2019, 08, 16)+ 1;  
                Veg_0 = 0; 
                Veg_max = 0.6; 
            elseif (20191007 <= str2num(Date)) && (str2num(Date) <= 20191231)  % 小麦
                Time_init = datenum(2019, 10, 07); 
                Total_days = datenum(2019, 12, 31) - datenum(2019, 10, 07)+ 1;  
                Veg_0 = 0; 
                Veg_max = 0.1; 
            else
                Time_init = 0;
                Total_days = 0;  
                Veg_0 = 0; 
                Veg_max = 0;
            end
        end
    elseif strcmp(Name_Tower, 'HL')   
        if str2num(Year) == 2017
            if (20170506 <= str2num(Date)) && (str2num(Date) <= 20171027)  % 玉米
                Time_init = datenum(2017, 05, 06); 
                Total_days = datenum(2017, 10, 27) - datenum(2017, 05, 06)+ 1;  
                Veg_0 = 0; 
                Veg_max = 2.5; 
            else
                Time_init = 0;
                Total_days = 0;
                Veg_0 = 0; 
                Veg_max = 0;
            end
        elseif str2num(Year) == 2018
            if (20180502 <= str2num(Date)) && (str2num(Date) <= 20181023)  % 玉米
                Time_init = datenum(2018, 05, 02); 
                Total_days = datenum(2018, 10, 23) - datenum(2018, 05, 02)+ 1;  
                Veg_0 = 0; 
                Veg_max = 2.5; 
            else
                Time_init = 0;
                Total_days = 0;  
                Veg_0 = 0; 
                Veg_max = 0;
            end
        elseif str2num(Year) == 2019                                       
            if (20190415 <= str2num(Date)) && (str2num(Date) <= 20190927)  % 杂草
                Time_init = datenum(2019, 04, 15); 
                Total_days = datenum(2019, 09, 27) - datenum(2019, 04, 15)+ 1;  
                Veg_0 = 0; 
                Veg_max = 1.0; 
            else
                Time_init = 0;
                Total_days = 0;  
                Veg_0 = 0; 
                Veg_max = 0;
            end
        end
    elseif strcmp(Name_Tower, 'DM')   
        if str2num(Year) == 2017
            if (20170513 <= str2num(Date)) && (str2num(Date) <= 20171101)  % 玉米
                Time_init = datenum(2017, 05, 13); 
                Total_days = datenum(2017, 11, 01) - datenum(2017, 05, 13)+ 1;  
                Veg_0 = 0; 
                Veg_max = 2.5; 
            else
                Time_init = 0;
                Total_days = 0;  
                Veg_0 = 0; 
                Veg_max = 0;
            end
        elseif str2num(Year) == 2018
            if (20180507 <= str2num(Date)) && (str2num(Date) <= 20181101)  % 玉米
                Time_init = datenum(2018, 05, 07); 
                Total_days = datenum(2018, 11, 01) - datenum(2018, 05, 07)+ 1;  
                Veg_0 = 0; 
                Veg_max = 2.5; 
            else
                Time_init = 0;
                Total_days = 0;  
                Veg_0 = 0; 
                Veg_max = 0;
            end
        elseif str2num(Year) == 2019                                       
            if (20190420 <= str2num(Date)) && (str2num(Date) <= 20191023)  % 玉米
                Time_init = datenum(2019, 04, 20); 
                Total_days = datenum(2019, 10, 23) - datenum(2019, 04, 20)+ 1;  
                Veg_0 = 0; 
                Veg_max = 2.5; 
            else
                Time_init = 0;
                Total_days = 0;  
                Veg_0 = 0; 
                Veg_max = 0;
            end
        end
    end
end