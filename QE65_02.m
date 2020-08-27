% The liuly Group 
% Key laboratory of Digital Earth Science, Aerospace Information Research
% Institute, Chinese Academy of Sciences, Beijing 100094, China
% 
% This example code illustrates how to access the observed spectral data 
% by QE65 spectrometer and the SIF retrievals
%
% If you have any questions, suggestions, comments on this example, please
% connect us (liuxj@radi.ac.cn; chenjidai@aircas.ac.cn)
% 
% Usage: run this script
%
% For more information, please read the program documentation! Please 
% refer to the relevant literature!
%
% Tested under: MATLAB R2019a
% Last updated: 2020-08-12

clc;                            % Clear the contents of the command window
clear;                          % Clear all variables from the workspace

% ----------主要功能: 大气校正---------------------------------------------------------- %%
% --此部分，注意修改
    Name_Tower = 'DM';                                    % 必修改: station
    Year = '2018';                                        % 必修改: year
% --路径设置
    Path_raw = strcat('.\Datasets\', Name_Tower, '\处理结果\', Year, '\', Name_Tower, '_Uncor_Cor_Rad_Ref\');       % 可修改：原始数据路径
    Path_pro = strcat('.\Datasets\', Name_Tower, '\处理结果\');                                                 % 可修改：处理结果保存路径
    Path_Meteo =  strcat('.\Settings\', Name_Tower, '_', Year, '_气压_温度_数据.txt');                      % 可修改：气压_温度数据路径
    
switch Name_Tower
    case 'XTS'
        % --设置上下行透过率路径
        Path_LUT_Dir = '.\Settings\Atmospheric_Transmittance\LLY_Ground_0050_SR_030_LUT_Dir_Tower_2.mat';      
        Path_LUT_Tot = '.\Settings\Atmospheric_Transmittance\LLY_Ground_0050_SR_030_LUT_Tot_Tower_2.mat';     
        % --对应站点的经纬度
        lon = 116.443;     lat = 40.179;
        TimeZone = 8;                     % 时区都为8,记录的北京时间
        % 海拔
        altitude = 50;
        % 塔高
        H_tower = 3.20;
        % 光纤观测地物角度
        Angle = 25;
        % 标准气压、温度
        P_zero = 1007.137;
        T_zero = 293.980;
    case 'HL'
        % --设置上下行透过率路径
        Path_LUT_Dir = '.\Settings\Atmospheric_Transmittance\LLY_Ground_0500_SR_030_LUT_Dir_Tower_2.mat';      
        Path_LUT_Tot = '.\Settings\Atmospheric_Transmittance\LLY_Ground_0500_SR_030_LUT_Tot_Tower_2.mat';    
        % --对应站点的经纬度
        lon = 115.788;     lat = 40.349;
        TimeZone = 8;
        % 海拔
        altitude = 500;
        % 塔高
        H_tower = 4.08;
        % 光纤观测地物角度
        Angle = 25;
        % 标准气压、温度
        P_zero = 955.887;
        T_zero = 293.980;
    case 'DM'   
        % --设置上下行透过率路径
        Path_LUT_Dir = '.\Settings\Atmospheric_Transmittance\LLY_Ground_1500_SR_030_LUT_Dir_Tower_20.mat';      
        Path_LUT_Tot = '.\Settings\Atmospheric_Transmittance\LLY_Ground_1500_SR_030_LUT_Tot_Tower_20.mat';    
        % --对应站点的经纬度
        lon = 100.372;     lat = 38.856;
        TimeZone = 8;
        % 海拔
        altitude = 1500;
        % 塔高
        H_tower = 24.5;
        % 光纤观测地物角度
        Angle = 25;
        % 标准气压、温度
        P_zero = 850.338;
        T_zero = 287.540;
        
    otherwise
        warning('Unexpected folder!')                                        % 报错：不存在此站点数据
end

%% ------读取未校正辐亮度数据、以及气压温度和透过率数据------ %%
% --加载大气校正相关数据：获取气压、温度数据以及上下行透过率
[Meteo] = Load_Meteo(Path_Meteo);                              % Load_Meteo是函数，加载气压、温度数据，为了计算等效长度
[Table_up, Table_down] = Load_LUT(Path_LUT_Dir, Path_LUT_Tot); % Load_LUT是函数，加载查找表；
% --加载未校正光谱数据，并进行大气校正
files = dir(strcat(Path_raw, Name_Tower,'_Uncor_', '*.xlsx'));  % 这里只读取了xlsx文件
nfiles = size(files, 1);                              % xlsx文件数量


for i = 1:1:nfiles
    Path_i = strcat(Path_raw, files(i).name);       % 文件路径
    [~, name, ext] = fileparts(files(i).name);        % 获取指定文件的路径、文件名和扩展名（文件格式）。
    Date = name(isstrprop(name, 'digit'));            % 获取字符串中的数字，即观测日期；strrep也可以实现  
    
    data = [];
    data_raw = xlsread(Path_i, 2);      
    num_meas = (size(data_raw, 2)-1)/2;
    num_wl = size(data_raw, 1)-2;
    data.wl = data_raw(3:num_wl+2,1);
    data.time = data_raw(1, 2:num_meas+1);
    data.sza = data_raw(2, 2:num_meas+1);
    data.veg = data_raw(3:num_wl+2, 2:num_meas+1);
    data.sky = data_raw(3:num_wl+2, num_meas+2:2*num_meas+1);
    data.nmeas = num_meas;
    data.npixels = num_wl;
    
    % --大气校正部分
    disp('Starting： 大气校正' + string(Date) + ', 之后保存');   tic;
    % --加载植被参数，高度，生长期，
    [Time_init, Total_days, Veg_0, Veg_max] = Setup_Veg(Name_Tower, Year, Date);      % 注意：相应站点的植被，观测天数，生长高度等
    veg_cor = [];
    sky_cor = [];
    for n = 1:1:data.nmeas
        % --植被的高度____注意：作物种植时间与作物最终高度 
        Year_i = str2num(strcat(Date(1), Date(2), Date(3), Date(4)));
        Month_i = str2num(strcat(Date(5), Date(6)));
        Day_i = str2num(strcat(Date(7), Date(8)));
        if Total_days == 0                              % 这里的判断，能够区分开
            H_veg = 0;
        else
            % 计算植被某个时间点的高度
            H_veg = ((datenum(Year_i, Month_i, Day_i) - Time_init) * (Veg_max-Veg_0) + Veg_0) ./ Total_days; 
        end
        
        % --上下行路径长度， 观测角度Angle及太阳天顶角SZA
        L_up = ((H_tower - H_veg) ./ 1000) ./ cosd(Angle);                 
        L_down = ((H_tower - H_veg) ./ 1000) ./ cosd(data.sza(1, n));
        % --等效路径
        Index_meteo = find(Meteo(:,1) == str2double(Date));
        length_meteo = size(Index_meteo, 1);
        if length_meteo == 0
            RTPL_up =  real(L_up * (P_zero./P_zero).^0.9353 * (T_zero./T_zero).^0.1936);
            RTPL_down =  real(L_down * (P_zero./P_zero).^0.9353 * (T_zero./T_zero).^0.1936);
        elseif length_meteo ~= 0
            Tem_interp = interp1(Meteo(Index_meteo(1,1):Index_meteo(end,1),2), Meteo(Index_meteo(1,1):Index_meteo(end,1),3), data.time(1,n), 'nearest');
            Pre_interp = interp1(Meteo(Index_meteo(1,1):Index_meteo(end,1),2), Meteo(Index_meteo(1,1):Index_meteo(end,1),4), data.time(1,n), 'nearest');
            RTPL_up =  real(L_up * (Pre_interp./P_zero).^0.9353 * ((T_zero-273.15)./Tem_interp).^0.1936);
            RTPL_down =  real(L_down * (Pre_interp./P_zero).^0.9353 * ((T_zero-273.15)./Tem_interp).^0.1936);
        end
        % --E790/E660, 660、790分别表示：660.044nm与790.113nm
        Index790 = data.wl == 790.113; 
        Index660 = data.wl == 660.044; 
        E_Ratio = data.sky(Index790, n) ./ data.sky(Index660, n);
        % --查找上下行的透过率, LUT_Tra函数
        [Up_tra] = LUT_Tra(Table_up, E_Ratio, RTPL_up);               % LUT_Tra是函数，查找透过率
        [Down_tra] = LUT_Tra(Table_down, E_Ratio, RTPL_down);
        New_Up_tra = [Up_tra(1,1:20), Up_tra,Up_tra(1,(end-20+1):end)];
        New_Down_tra = [Down_tra(1,1:20), Down_tra, Down_tra(1,(end-20+1):end)];
        [~, Up_tra_min_num] = min(New_Up_tra(1,:), [], 2);
        [~, ~] = min(New_Down_tra(1,:), [], 2);
        Num_tra = Up_tra_min_num;                                          % 上行透过率的最小值更稳定
        Window_O2A = [758, 770]; 
        [Range, Index_veg_min, Index_sky_min] = Cal_Index(Window_O2A, data); % Cal_Index为函数，找到地物和入照曲线，氧气A波段最小值所在序号
        Num_qe = (Index_veg_min + Index_sky_min) ./ 2;
        Inter_bef = Num_tra - (Num_qe - 1);
        Inter_aft = Num_tra + (1021 - Num_qe);
        
        % --将TOA(传感器处)数值转化为TOC(冠层处)数值-----
        veg_cor(:, n) = data.veg(:,n) ./ New_Up_tra(1,Inter_bef:Inter_aft)';
        sky_cor(:, n) = data.sky(:,n) .* New_Down_tra(1,Inter_bef:Inter_aft)';
    end
    
    data.veg_cor = veg_cor;
    data.sky_cor = sky_cor;
    data.ref_cor = data.veg_cor ./ data.sky_cor;

    % --保存数据为nc, xlsx和mat格式
    save_xlsx(Path_pro, Name_Tower, Date, data);
    save_mat(Path_pro, Name_Tower, Date, data);
    save_nc(Path_pro, Name_Tower, Date, data);
    
end
disp('Time delays: ' + string(toc) + 's');


%% ------------------------函数部分------------------------ %%
% --函数备注：判断窗口内入照和地物光谱最小值和序列号
function [Range, Index_veg_min_ave, Index_sky_min_ave] = Cal_Index(wl_range, data)
    Range = find(data.wl > wl_range(1) & data.wl < wl_range(2));    % 吸收波段范围的索引
    Index_left = Range(1);
    Index_right = Range(end);
    for i = 1:1:data.nmeas
        [veg_min(1,i), Index_veg_min(1,i)] = min(data.veg(Index_left:Index_right, i),[],1);
        [sky_min(1,i), Index_sky_min(1,i)] = min(data.sky(Index_left:Index_right, i),[],1);
    end
    % 入照、地物的最小值在初晨、傍晚的数据存在异常，除去初晨与傍晚各10%有问题的数据 %
    Bound_min = round(data.nmeas*0.1);
    Bound_max = round(data.nmeas)-round(data.nmeas*0.1);
    Index_veg_min_ave = round(mean(Index_veg_min(1,Bound_min:Bound_max),2))+ Index_left-1;
    Index_sky_min_ave = round(mean(Index_sky_min(1,Bound_min:Bound_max),2))+ Index_left-1;
end

% --函数备注：加载气压、温度数据，并格式转换。如：2017/5/1 0:20转换为：20170501 0.01389；
function [Meteo] = Load_Meteo(Path_Meteo)

    [Date, Time, Tem, Press] = textread(Path_Meteo, '%s%s%s%s', 'headerlines', 1);
    % ----------“年份+日期”格式转换---------- %
    nyear = size(Date, 1);
    for i = 1:1:nyear
        Date_i = Date{i,1};
        Length_Date_i = length(Date_i);
        Num_dia = strfind(Date_i,'/');
        % 提取年数据
        Year = strcat(Date_i(1), Date_i(2), Date_i(3), Date_i(4));
        % 提取月数据
        if Num_dia(1,2) - Num_dia(1,1) == 2
            Month = strcat('0', Date_i(Num_dia(1,1)+1));
        elseif Num_dia(1,2) - Num_dia(1,1) == 3
            Month = strcat(Date_i(Num_dia(1,1)+1), Date_i(Num_dia(1,1)+2));
        end
        % 提取日数据
        if Length_Date_i - Num_dia(1,2) == 1
            Day = strcat('0', Date_i(Num_dia(1,2)+1));
        elseif Length_Date_i - Num_dia(1,2) == 2
            Day = strcat(Date_i(Num_dia(1,2)+1), Date_i(Num_dia(1,2)+2));
        end
        Date_new(i,1) = str2double(strcat(Year, Month, Day));
    end
    
    % ----------“时+分”格式转换---------- %
    ntime = size(Time, 1);
    for j = 1:1:ntime
        Time_j = Time{j, 1};
        Num_col = strfind(Time_j, ':');
        %提取时、分数据
        if Num_col(1,1) == 2
            Hour = str2double(Time_j(Num_col(1,1)-1));
            Minute = str2double(strcat(Time_j(Num_col(1,1)+1), Time_j(Num_col(1,1)+2)));
        elseif Num_col(1,1) == 3
            Hour = str2double(strcat(Time_j(Num_col(1,1)-2), Time_j(Num_col(1,1)-1)));
            Minute = str2double(strcat(Time_j(Num_col(1,1)+1), Time_j(Num_col(1,1)+2)));
        end
        Time_new(j,1) = Hour./24 + Minute./60./24;
    end
    
    % ----------温度、气压：字符格式转换为数字格式---------- %
    ntem = size(Tem, 1);
    for k = 1:1:ntem
        Tem_new(k, 1) = str2double(Tem{k,1});
        Press_new(k, 1) = str2double(Press{k,1});
    end

    Meteo = [Date_new, Time_new, Tem_new, Press_new];
end

% --函数备注：加载上下行透过率，MORTRAN模拟数据
function [Table_up_new, Table_down_new] = Load_LUT(Path_LUT_Dir, Path_LUT_Tot)
    disp('Starting： 加载上下行透过率 ');
    % ----读取E比值 _Tdir的模拟数据 && _Ttot的模拟数据---- %
    LUT_Dir = load(Path_LUT_Dir);      Name_Dir = whos('-file', Path_LUT_Dir);  % 不用知道结构名来访问数据
    LUT_Tot = load(Path_LUT_Tot);      Name_Tot = whos('-file', Path_LUT_Tot); 
    Table_up_raw = LUT_Dir.(Name_Dir.name);
    Table_down_raw = LUT_Tot.(Name_Tot.name); 
    AOD_00_90_up = Table_up_raw(:,2);  
    AOD_00_90_down = Table_down_raw(:,2);
    for i = 1:1:10             % 对上行\下行 模拟数据按照AOD固定，VZA由小到大排序
        AOD_00_90_Num_up(:,i) = find(AOD_00_90_up(:,1) == i./10);
        AOD_00_90_Num_down(:,i) = find(AOD_00_90_down(:,1) == i./10);
    end
    Table_up_new =[];
    Table_down_new =[];
    for j = 1:1:10
        for k = 1:1:180
            Table_up_new = [Table_up_new; Table_up_raw(AOD_00_90_Num_up(k,j), :)];
            Table_down_new = [Table_down_new; Table_down_raw(AOD_00_90_Num_down(k,j), :)];
        end
    end
    
end

% ---函数备注：通过查找表的方式，计算得到上行透过率
function [Transmittance] = LUT_Tra(Table_up_new, E_Ratio, RTPL_up)

    % ---当H_tower/cos(VZA)确定，对E比值、Tdir透过率的数值进行线性插值
    L_min = min(Table_up_new(:, 1));
    L_max = max(Table_up_new(:, 1));
    if (RTPL_up >= L_min && RTPL_up <= L_max)
        for i = 1:1:10
            Up_value_interp(i, 1:1023) = interp1(Table_up_new((180*i-179):180*i,1), Table_up_new((180*i-179):180*i,4:1026), RTPL_up, 'linear');
        end
    elseif (RTPL_up < L_min)
        RTPL_up_min = L_min;
        for i = 1:1:10
            Up_value_interp(i, 1:1023) = interp1(Table_up_new((180*i-179):180*i,1), Table_up_new((180*i-179):180*i,4:1026), RTPL_up_min, 'linear');
        end
    elseif (RTPL_up > L_max)
        RTPL_up_max = L_max;
        for i = 1:1:10
            Up_value_interp(i, 1:1023) = interp1(Table_up_new((180*i-179):180*i,1), Table_up_new((180*i-179):180*i,4:1026), RTPL_up_max, 'linear');
        end
    end

    % ---当E比值确定，对Tdir透过率的数值进行线性插值
    E_min = min(Up_value_interp(:, 1));
    E_max = max(Up_value_interp(:, 1));
    if (E_Ratio >= E_min && E_Ratio <= E_max)
        Up_tra_interp(1,:) = interp1(Up_value_interp(1:10,1), Up_value_interp(1:10,3:1023), E_Ratio, 'linear');
    elseif (E_Ratio < E_min)
        E_Ratio_min = E_min;
        Up_tra_interp(1,:) = interp1(Up_value_interp(1:10,1), Up_value_interp(1:10,3:1023), E_Ratio_min, 'linear');
    elseif (E_Ratio > E_max)
        E_Ratio_max = E_max;
        Up_tra_interp(1,:) = interp1(Up_value_interp(1:10,1), Up_value_interp(1:10,3:1023), E_Ratio_max, 'linear');
    end
    Transmittance(1, :) = Up_tra_interp;
end

% ---函数备注：保存数据
function save_xlsx(Path_pro, Name_Tower, Date, data)
    disp('Saving： xlsx');  
    
    % --------------建立文件夹-------------- %
    folder_date = strcat(Path_pro, Date(1:4), '\');    
    if exist(folder_date, 'dir')==0     % 判断【年份】文件夹是否存在
        mkdir(folder_date);             % 不存在时候，创建文件夹
    end
    folder = strcat(folder_date, Name_Tower, '_Uncor_Cor_Rad_Ref\');     
    if exist(folder, 'dir')==0          % 判断【_Uncor_Rad_Ref_】文件夹是否存在
        mkdir(folder);                  % 不存在时候，创建文件夹
    end
    % --------------输出Excel结果-------------- %
    Output_path = strcat(folder, Name_Tower, '_Cor_Rad_Ref_', Date, '.xlsx');       % 修改:输出.xlsx文件夹的观测站名称，如：DM
    Out_Ref_time_sza = [data.time; data.sza];
    Out_Rad_time_sza = [data.time, data.time; data.sza, data.sza];
    Out_rad_veg_sky = [data.veg_cor, data.sky_cor];
    Names = {'时间'; 'SZA'};
    % --地物与入照的辐亮度（未校正）-- %
    xlswrite(Output_path, Names, '地物与入照的辐亮度（已校正）', 'A1:A2');
    xlswrite(Output_path, data.wl, '地物与入照的辐亮度（已校正）', 'A3');
    xlswrite(Output_path, Out_Rad_time_sza, '地物与入照的辐亮度（已校正）', 'B1');
    xlswrite(Output_path, Out_rad_veg_sky, '地物与入照的辐亮度（已校正）', 'B3');
    % --反射率（未校正）-- %
    xlswrite(Output_path, Names, '反射率（已校正）', 'A1:A2');
    xlswrite(Output_path, data.wl, '反射率（已校正）', 'A3');
    xlswrite(Output_path, Out_Ref_time_sza, '反射率（已校正）', 'B1');
    xlswrite(Output_path, data.ref_cor, '反射率（已校正）', 'B3');
    
end
function save_mat(Path_pro, Name_Tower, Date, data)
    disp('Saving： mat'); 

    % --------------建立文件夹-------------- %
    folder_date = strcat(Path_pro, Date(1:4), '\');    
    if exist(folder_date, 'dir')==0     
        mkdir(folder_date);            
    end
    folder = strcat(folder_date, Name_Tower, '_Uncor_Cor_Rad_Ref\');  
    if exist(folder, 'dir')==0    
        mkdir(folder);             
    end
    % --------------输出mat结果-------------- %
    Output_path = strcat(folder, Name_Tower, '_Cor_Rad_Ref_', Date, '.mat');       % 输出
    
    wl = data.wl;
    time = data.time;
    sza = data.sza;
    veg = data.veg_cor;
    sky = data.sky_cor;
    ref = data.ref_cor;
    save(Output_path, 'wl', 'time', 'sza', 'veg', 'sky', 'ref');
end
function save_nc(Path_pro, Name_Tower, Date, data)
    disp('Saving： nc'); 

    % --------------建立文件夹-------------- %
    folder_date = strcat(Path_pro, Date(1:4), '\');    
    if exist(folder_date, 'dir')==0     % 判断【年份】文件夹是否存在
        mkdir(folder_date);             % 不存在时候，创建文件夹
    end
    folder = strcat(folder_date, Name_Tower, '_Uncor_Cor_Rad_Ref\');     
    if exist(folder, 'dir')==0          % 判断【_Uncor_Rad_Ref_】文件夹是否存在
        mkdir(folder);                  % 不存在时候，创建文件夹
    end
    % --------------输出nc结果-------------- %
    Output_path = strcat(folder, Name_Tower, '_Cor_Rad_Ref_', Date, '.nc');       % 修改:输出.xlsx文件夹的观测站名称，如：DM
    
    % ----------DEFINE THE FILE---------- %                                        
    ncid = netcdf.create(Output_path,'CLOBBER');   % 创建一个存放数据的nc文件
    %-----------define dimension--------- %   
    dimidx = netcdf.defDim(ncid,'Wavelength',data.npixels);
    dimidy = netcdf.defDim(ncid,'Measures', data.nmeas);    
    
    %-----------define new variables------ %
    varid1 = netcdf.defVar(ncid,'Time','double',[dimidy]);
    varid2 = netcdf.defVar(ncid,'SZA','double',[dimidy]);
    varid3 = netcdf.defVar(ncid,'wl','double',[dimidx]);
    varid4 = netcdf.defVar(ncid,'veg','double',[dimidx dimidy]);
    varid5 = netcdf.defVar(ncid,'sky','double',[dimidx dimidy]);
    varid6 = netcdf.defVar(ncid,'ref','double',[dimidx dimidy]);
    
    %------------define attributes of the new variables---------%  
    netcdf.putAtt(ncid,varid1,'units','/');       %单位信息和long_name，其它的信息可依此定义
    netcdf.putAtt(ncid,varid2,'units','degree');      
    netcdf.putAtt(ncid,varid3,'units','nm');  
    netcdf.putAtt(ncid,varid4,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid5,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid6,'units','mW/m2/nm/sr');
   
    netcdf.putAtt(ncid,varid1,'long_name','decimal Julian day');        
    netcdf.putAtt(ncid,varid2,'long_name','Solar zenith angle');   
    netcdf.putAtt(ncid,varid3,'long_name','Wavelength');  
    netcdf.putAtt(ncid,varid4,'long_name','Corrected Canopy radiance');
    netcdf.putAtt(ncid,varid5,'long_name','Corrected Solar irradiance / pi');
    netcdf.putAtt(ncid,varid6,'long_name','Corrected Canopy reflectance');

    netcdf.endDef(ncid);
    % 放入变量
    netcdf.putVar(ncid,varid1,data.time);
    netcdf.putVar(ncid,varid2,data.sza);
    netcdf.putVar(ncid,varid3,data.wl);
    netcdf.putVar(ncid,varid4,data.veg_cor);
    netcdf.putVar(ncid,varid5,data.sky_cor);
    netcdf.putVar(ncid,varid6,data.ref_cor);
    
    netcdf.close(ncid);
end
