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

% ----------主要功能: organize the observed radiance data in txt format into xlsx and mat format---------- %%
% --此部分，注意修改
    Name_Tower = 'DM';                                    % 必修改: station
    Year = '2018';                                        % 必修改: year
% --路径设置
    Path_raw = strcat('.\Datasets\', Name_Tower, '\原始数据\', Year, '\');    % 可修改：原始数据路径（事先建立此文件夹, 把原始数据放在其中）
    Path_pro = strcat('.\Datasets\', Name_Tower, '\处理结果\');               % 可修改：处理结果保存路径（事先建立此文件夹）
    Norm_wl(:, 1) = importdata('.\Settings\wl_pro.txt');                 % 可修改：标准波长路径
    
switch Name_Tower
    case 'XTS'
        % --对应站点的仪器序号
        QEPB_serial = 'QEPB1220';
        if strcmp(Name_Tower, 'XTS') && strcmp(Year,'2017')                      %  --***注意：2017年小汤山和怀来仪器,7月6号---- %                                             
            QEPB_serial = 'QEPB1218';
        end
        % --对应站点的经纬度
        lon = 116.443;     lat = 40.179;
        TimeZone = 8;                     % 时区都为8,记录的北京时间
    case 'HL'     
        % --对应的仪器序号
        QEPB_serial = 'QEPB1218';
        if strcmp(Name_Tower, 'HL') && strcmp(Year,'2017')                        %  --***注意：2017年小汤山和怀来仪器---- %                                             
            QEPB_serial = 'QEPB1220';
        elseif strcmp(Name_Tower, 'HL') && strcmp(Year,'2020')                    %  --***注意：2020年阿柔和怀来仪器更换---- %                                             
            QEPB_serial = 'QEPB1564';
        end
        % --对应站点的经纬度
        lon = 115.788;     lat = 40.349;
        TimeZone = 8;
    case 'DM'            
        % --对应的仪器序号
        QEPB_serial = 'QEPB1219';
        % --对应站点的经纬度
        lon = 100.372;     lat = 38.856;
        TimeZone = 8;

    otherwise
        warning('Unexpected folder!')                                        % 报错：不存在此站点数据
end

%% ------加载观测光谱_txt格式文本文件，以标准波段插值------ %%
Files = dir(fullfile(Path_raw, strcat(Year, '*')));  % Struct of files
Length = size(Files, 1);                             % Number of daily observation data,days

for i = 1:1:Length
    
    date_i = Files(i,1).name;
    year_i = str2double(strcat(date_i(1), date_i(2), date_i(3), date_i(4)));
    month_i = str2double(strcat(date_i(5), date_i(6)));
    day_i = str2double(strcat(date_i(7), date_i(8)));
    if strcmp(Name_Tower, 'XTS') && (str2num(date_i) >= 20170706)  % 更换仪器
        QEPB_serial = 'QEPB1220';
    elseif strcmp(Name_Tower, 'HL') && (str2num(date_i) >= 20170711)  % 更换仪器
        QEPB_serial = 'QEPB1218';
    end
    path_i = strcat(Path_raw, Files(i,1).name, '\Auto\', QEPB_serial, '\'); % 最终文件夹路径
    disp('Starting： 读取' + string(date_i) + '日期所有txt文件，并以标准波段插值, 之后保存');   tic;
    %--读取最终文件夹的txt数据
    txts = dir(fullfile(path_i, '*.txt'));     % 当天所有的txt文件，结构体形式
    ntxts = size(txts, 1);                     % 当天，txt的数量   
    wl = [];                                   % **----这里很重要，清零数据，不然会有内存
    time = [];
    sza = [];
    spec = [];
    for j = 1:1:ntxts
        path_j = strcat(path_i, txts(j, 1).name);                           % 最终txt格式文件路径
        [wl_j, spec_j, time_j] = get_spec(path_j);
        sza_j = Cal_SZA(time_j, year_i, month_i, day_i, lon, lat, TimeZone);
        wl(:,j) = wl_j(:,1);
        sza(:,j) = sza_j(:,1);
        time(:,j) = time_j(:,1);
        spec(:,j) = spec_j(:,1);
    end
    % --对数据按照时间升序排序-------- %
    [time, index] = sort(time);
    sza = sza(index);
    spec = spec(:, index);
    % --对地物、入照数据整理-------- %
    time_mean = [];
    sza_mean = [];
    veg = [];
    sky = [];
    % --判断txt数据的数量是否为3的倍数-------- %
    if mod(ntxts, 3) == 0                   
        for k = 1:1:(ntxts/3)
            time_mean(:, k) = mean(time(1, (3*k-2):3*k), 2);
            sza_mean(:, k) = mean(sza(1, (3*k-2):3*k), 2);
            veg(:, k) = spec(:, 3*k-1);
            % --剔除入照有问题的数据-- %,    --{质量控制内容（后期更新}--
            % ***注意修改***，这里好像效果不好【Irr_760/Irr_758】
            S1 = mean(spec(674:676, 3*k-2), 1) ./ mean(spec(657:659, 3*k-2), 1);
            S3 = mean(spec(674:676, 3*k), 1) ./ mean(spec(657:659, 3*k), 1);
            if  le(S1, 0.5) && le(S3, 0.5)
                sky(:,k) = (spec(:,3*k-2) + spec(:,3*k))./2.0;
            elseif gt(S1, 0.5) && gt(S3, 0.5)
                sky(:,k) = (spec(:,3*k-2) + spec(:,3*k))./2.0;
            elseif gt(S1, 0.5) && le(S3, 0.5)
                sky(:,k) = spec(:,3*k);
            elseif le(S1, 0.5) && gt(S3, 0.5)
                sky(:,k) = spec(:,3*k-2);
            end
        end
    else
        warning = 'The number of txt is wrong!';
    end
    % ---以标准波段插值---------------%
    Prim_wl(:,1) = wl(:,1);
    veg_interp = [];
    sky_interp = [];
    for m = 1:1:(ntxts/3)  
        veg_interp(:,m) = interp1(Prim_wl(:, 1), veg(:, m), Norm_wl(:,1), 'spline');
        sky_interp(:,m) = interp1(Prim_wl(:, 1), sky(:, m), Norm_wl(:,1), 'spline');
    end
    % **注意数据___【张掖：20170511-20170711为特殊数据，需要转换单位，即入照、地物数值单位转换到太阳辐亮度、地物辐亮度】----%
    if (str2double(date_i) >= 20170511) && (str2double(date_i) <= 20170711)
        veg_interp = veg_interp*10;
        sky_interp = sky_interp*10/pi;    % 注意***，这里定标系数有问题（by胡娇婵）
    end
    
    data = [];
    data.wl = Norm_wl;
    data.time = time_mean;
    data.npixels = size(data.wl, 1);
    data.nmeas = size(data.time, 2);
    data.sza = sza_mean;
    data.veg = veg_interp;
    data.sky = sky_interp;
    data.ref = data.veg ./ data.sky;
    
    % --保存数据为nc, xlsx和mat格式
    save_xlsx(Path_pro, Name_Tower, date_i, data);
    save_mat(Path_pro, Name_Tower, date_i, data);
    save_nc(Path_pro, Name_Tower, date_i, data);
    
    disp('Time delays: ' + string(toc) + 's');
end

%% ------------------------函数部分------------------------ %%
% --函数备注：读取光谱数据
function [wl, spec, time] = get_spec(path)
    % --读取时间
    if nargin == 1                                         % 判断输入变量的个数：
        fid = fopen(path);                                 % fid为fopen命令返回的文件标识符
        a = textscan(fid, '%s%s%s%s', 'headerlines',10);   % 找到starting
        fclose(fid);
        a4 = a{1, 4};
        T = a4{1, 1};
        L = length(T);
        
        if L == 7
            Hour = str2double(T(1));
            Min = str2double(strcat(T(3), T(4)));
            Sec = str2double(strcat(T(6), T(7)));
        elseif L == 8
            Hour = str2double(strcat(T(1), T(2)));
            Min = str2double(strcat(T(4), T(5)));
            Sec = str2double(strcat(T(7), T(8)));
        end
        time = (Hour/24.0 + Min/24.0/60.0 + Sec/24.0/60.0/60.0);   % 将时间换算成0-1的时间范围
    end
    
    % --读取数据
    fid = fopen(path);
    b = textscan(fid, '%s%s%s%s', 'headerlines', 18);
    fclose(fid);
    b1 = b{1, 1};
    b4 = b{1, 4};
    N = size(b1,1);
    
    wl = str2num(char(b1));    % str2num,而str2double只适用于标量的转换spec=str2double(b1);
    spec = str2num(char(b4));     
end

% --函数备注：保存数据为xlsx和mat
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
    Output_path = strcat(folder, Name_Tower, '_Uncor_Rad_Ref_', Date, '.xlsx');       % 修改:输出.xlsx文件夹的观测站名称，如：DM
    
    Out_Ref_time_sza = [data.time; data.sza];
    Out_Rad_time_sza = [data.time, data.time; data.sza, data.sza];
    Out_rad_veg_sky  = [data.veg, data.sky];
    Names = {'时间'; 'SZA'};
    % --地物与入照的辐亮度（未校正）-- %
    xlswrite(Output_path, Names, '地物与入照的辐亮度（未校正）', 'A1:A2');
    xlswrite(Output_path, data.wl, '地物与入照的辐亮度（未校正）', 'A3');
    xlswrite(Output_path, Out_Rad_time_sza, '地物与入照的辐亮度（未校正）', 'B1');
    xlswrite(Output_path, Out_rad_veg_sky, '地物与入照的辐亮度（未校正）', 'B3');

    % --反射率（未校正）-- %
    xlswrite(Output_path, Names, '反射率（未校正）', 'A1:A2');
    xlswrite(Output_path, data.wl, '反射率（未校正）', 'A3');
    xlswrite(Output_path, Out_Ref_time_sza, '反射率（未校正）', 'B1');
    xlswrite(Output_path, data.ref, '反射率（未校正）', 'B3');

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
    Output_path = strcat(folder, Name_Tower, '_Uncor_Rad_Ref_', Date, '.mat');       % 输出
    
    wl = data.wl;
    time = data.time;
    sza = data.sza;
    veg = data.veg;
    sky = data.sky;
    ref = data.ref;
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
    Output_path = strcat(folder, Name_Tower, '_Uncor_Rad_Ref_', Date, '.nc');       % 修改:输出.xlsx文件夹的观测站名称，如：DM
    
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
    netcdf.putAtt(ncid,varid4,'long_name','Canopy radiance');
    netcdf.putAtt(ncid,varid5,'long_name','Solar irradiance / pi');
    netcdf.putAtt(ncid,varid6,'long_name','Canopy reflectance');

    netcdf.endDef(ncid);
    % 放入变量
    netcdf.putVar(ncid,varid1,data.time);
    netcdf.putVar(ncid,varid2,data.sza);
    netcdf.putVar(ncid,varid3,data.wl);
    netcdf.putVar(ncid,varid4,data.veg);
    netcdf.putVar(ncid,varid5,data.sky);
    netcdf.putVar(ncid,varid6,data.ref);
    
    netcdf.close(ncid);
end

