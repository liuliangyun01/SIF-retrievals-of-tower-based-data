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

% ----------主要功能: 整合数据，并30分钟平均---------------------------------------------------------- %%
% --此部分，注意修改
    Name_Tower = 'HL';                                    % 必修改: station
    Year = '2018';                                        % 必修改: year
% --路径设置
    
    uncor_cor{1, 1} = char('Uncor');      % 对大气校正前后的数据，进行整理；
    uncor_cor{2, 1} = char('Cor'); 
    % 可修改：未校正辐亮度数据路径
    Path_Uncor_Cor = strcat('.\Datasets\', Name_Tower, '\处理结果\', Year, '\', Name_Tower, '_Uncor_Cor_SIF_Irr_Rad_Refl\');
    % 可修改：结果存放路径,未校正和校正
    Path_pro = strcat('.\Datasets\', Name_Tower, '\处理结果\', Year, '\', Name_Tower, '_Uncor_Cor_30min_or_not\');  
    Norm_wl(:, 1) = importdata('.\Settings\wl_pro.txt');                 % 可修改：标准波长路径
    
switch Name_Tower
    case 'XTS'
        % --对应站点的经纬度
        lon = 116.443;     lat = 40.179;
        TimeZone = 8;                     % 时区都为8,记录的北京时间
    case 'HL'     
        % --对应站点的经纬度
        lon = 115.788;     lat = 40.349;
        TimeZone = 8;
    case 'DM'            
        % --对应站点的经纬度
        lon = 100.372;     lat = 38.856;
        TimeZone = 8;
        
    otherwise
        warning('Unexpected folder!')                                        % 报错：不存在此站点数据
end


%% ------读取未校正辐亮度数据、以及气压温度和透过率数据------ %%
for i = 1:1:2
    name_uncor_cor = uncor_cor{i, 1};
    Path_raw = Path_Uncor_Cor;
    folder = Path_pro;           
    if exist(folder, 'dir')==0    
        mkdir(folder);            
    end
    Path_pro = folder; 
    
    % ----获取未校正和校正后的辐亮度数据，并大反演SIF---- %
    files = dir(strcat(Path_raw, Name_Tower, '_', name_uncor_cor, '*.nc'));   % 这里只读取了nc格式文件
    num_files = size(files, 1);                                               % nc文件数量
    allData=[];  allDate = []; allTime=[];  allSZA=[];
    aveData=[];  stdData=[];   numData=[];
    Date_30min=[];
    Time_30min=[];
    SZA_30min=[];
    for j = 1:1:num_files
        ncfilename=strcat(Path_raw, files(j).name);
        % 时间，sza, wavelength
        Time = ncread(ncfilename, 'Time');
        SZA = ncread(ncfilename, 'SZA');
        wl = ncread(ncfilename, 'wl');
        % 整理辐亮度数据
        Rad_680 = ncread(ncfilename, 'Rad_680'); Rad_705 = ncread(ncfilename, 'Rad_705'); Rad_740 = ncread(ncfilename, 'Rad_740'); Rad_770 = ncread(ncfilename, 'Rad_770'); Rad_800 = ncread(ncfilename, 'Rad_800');
        Irr_680 = ncread(ncfilename, 'Irr_680'); Irr_705 = ncread(ncfilename, 'Irr_705'); Irr_740 = ncread(ncfilename, 'Irr_740'); Irr_770 = ncread(ncfilename, 'Irr_770'); Irr_800 = ncread(ncfilename, 'Irr_800');
        Refl_680 = ncread(ncfilename, 'Refl_680');   Refl_705 = ncread(ncfilename, 'Refl_705');   Refl_740 = ncread(ncfilename, 'Refl_740');   Refl_770 = ncread(ncfilename, 'Refl_770');   Refl_800 = ncread(ncfilename, 'Refl_800');
        % 整理Ha
        Ha_sFLD = ncread(ncfilename, 'SIF_Ha_sFLD');
        Ha_3FLD = ncread(ncfilename, 'SIF_Ha_3FLD');
        Ha_iFLD = ncread(ncfilename, 'SIF_Ha_iFLD');
        Ha_pFLD = ncread(ncfilename, 'SIF_Ha_pFLD');
        Ha_SFM = ncread(ncfilename, 'SIF_Ha_SFM');
        Ha_SVD = ncread(ncfilename, 'SIF_Ha_SVD');
        Ha_DOAS = ncread(ncfilename, 'SIF_Ha_DOAS');
        % 整理O2B
        O2B_sFLD = ncread(ncfilename, 'SIF_O2B_sFLD');
        O2B_3FLD = ncread(ncfilename, 'SIF_O2B_3FLD');
        O2B_iFLD = ncread(ncfilename, 'SIF_O2B_iFLD');
        O2B_pFLD = ncread(ncfilename, 'SIF_O2B_pFLD');
        O2B_SFM = ncread(ncfilename, 'SIF_O2B_SFM');
        O2B_SVD = ncread(ncfilename, 'SIF_O2B_SVD');
        O2B_DOAS = ncread(ncfilename, 'SIF_O2B_DOAS');
        % 整理H2O
        H2O_sFLD = ncread(ncfilename, 'SIF_H2O_sFLD');
        H2O_3FLD = ncread(ncfilename, 'SIF_H2O_3FLD');
        H2O_iFLD = ncread(ncfilename, 'SIF_H2O_iFLD');
        H2O_pFLD = ncread(ncfilename, 'SIF_H2O_pFLD');
        H2O_SFM = ncread(ncfilename, 'SIF_H2O_SFM');
        H2O_SVD = ncread(ncfilename, 'SIF_H2O_SVD');
        H2O_DOAS = ncread(ncfilename, 'SIF_H2O_DOAS');
        % 整理O2A
        O2A_sFLD = ncread(ncfilename, 'SIF_O2A_sFLD');
        O2A_3FLD = ncread(ncfilename, 'SIF_O2A_3FLD');
        O2A_iFLD = ncread(ncfilename, 'SIF_O2A_iFLD');
        O2A_pFLD = ncread(ncfilename, 'SIF_O2A_pFLD');
        O2A_SFM = ncread(ncfilename, 'SIF_O2A_SFM');
        O2A_SVD = ncread(ncfilename, 'SIF_O2A_SVD');
        O2A_DOAS = ncread(ncfilename, 'SIF_O2A_DOAS');
        % 整理全波段
        Fullband = ncread(ncfilename, 'SIF_Fullband');
        
        DOY = floor(Time);
        Clock = Time - DOY;
        Date_j = files(j).name(isstrprop(files(j).name, 'digit'));
        Date = repmat(str2num(Date_j), size(Clock, 1), 1);     % 生成日期
        year = str2num(Date_j(1:4));
        month = str2num(Date_j(5:6));
        day = str2num(Date_j(7:8));  
        
        DataToAve = [Rad_680, Rad_705, Rad_740, Rad_770, Rad_800, ...
                     Irr_680, Irr_705, Irr_740, Irr_770, Irr_800, ...
                    Ha_sFLD,  Ha_3FLD,  Ha_iFLD,  Ha_pFLD,  Ha_SFM,  Ha_SVD,  Ha_DOAS, ...
                    O2B_sFLD, O2B_3FLD, O2B_iFLD, O2B_pFLD, O2B_SFM, O2B_SVD, O2B_DOAS, ...
                    H2O_sFLD, H2O_3FLD, H2O_iFLD, H2O_pFLD, H2O_SFM, H2O_SVD, H2O_DOAS, ...
                    O2A_sFLD, O2A_3FLD, O2A_iFLD, O2A_pFLD, O2A_SFM, O2A_SVD, O2A_DOAS, ...
                    Fullband];
        allData = [allData; DataToAve];
        allDate = [allDate; Date];
        allTime = [allTime; Time];
        allSZA = [allSZA; SZA];

        time_step=1/48;
        for k=1:48
            TimeWindow_LB = k*time_step-time_step;
            TimeWindow_UB = k*time_step;
            index_inTimeWindow = find((Clock>TimeWindow_LB & Clock<=TimeWindow_UB));
            DataToAve_Window = DataToAve(index_inTimeWindow, :);
            if size(DataToAve_Window,1)>0
                Q1 = prctile(DataToAve_Window,25,1);
                %Q2=prctile(DataToAve_Window,50,1);
                Q3 = prctile(DataToAve_Window,75,1);
                IQR = Q3-Q1;
                %根据Q1-1.5IQR, Q3+1.5IQR原则确定半小时内正常值上下边界（箱线图）
                lbound = Q1-1.5.*IQR;
                ubound = Q3+1.5.*IQR;

                for v=1:size(DataToAve_Window,2)     % v代表多少个类别
                    [index_valid,colom]=find(DataToAve_Window(:,v)>lbound(v) & DataToAve_Window(:,v)<ubound(v));
                    DataToAve_valid = DataToAve_Window(index_valid,v);
                    aveData_v(1,v) = mean(DataToAve_valid, 'omitnan');
                    stdData_v(1,v) = std(DataToAve_valid, 'omitnan');
                    numData_v(1,v) = length(DataToAve_valid);
                end
                aveData = [aveData; aveData_v];
                stdData = [stdData; stdData_v];
                numData = [numData; numData_v];
                Time_30min = [Time_30min; k*time_step+DOY(1)];
                SZA_30min_k = Cal_SZA(k*time_step, year, month, day, lon, lat, TimeZone);
                SZA_30min = [SZA_30min;SZA_30min_k];
                Date_30min_k = str2num(Date_j);     % 生成日期
                Date_30min = [Date_30min; Date_30min_k];
            end
        end 
    end

    % --保存数据为nc, xlsx和mat格式（不平均数据or30min平均数据）
    mean_or_not{1, 1} = char('all');      
    mean_or_not{2, 1} = char('30min');  
    for k = 1:1:2
        name_mean_or_not = mean_or_not{k, 1};
        if strcmp(name_mean_or_not, 'all')
%             save_xlsx(Path_pro, Name_Tower, name_uncor_cor, name_mean_or_not, Year, allData, allDate, allTime, allSZA, Norm_wl);
            save_mat(Path_pro, Name_Tower, name_uncor_cor, name_mean_or_not, Year, allData, allDate, allTime, allSZA, Norm_wl);
            save_nc(Path_pro, Name_Tower, name_uncor_cor, name_mean_or_not, Year, allData, allDate, allTime, allSZA, Norm_wl);
        elseif strcmp(name_mean_or_not, '30min')
%             save_xlsx(Path_pro, Name_Tower, name_uncor_cor, name_mean_or_not, Year, allData, allDate, allTime, allSZA, Norm_wl);
            save_mat(Path_pro, Name_Tower, name_uncor_cor, name_mean_or_not, Year, aveData, Date_30min, Time_30min, SZA_30min, Norm_wl);
            save_nc(Path_pro, Name_Tower, name_uncor_cor, name_mean_or_not, Year, aveData, Date_30min, Time_30min, SZA_30min, Norm_wl);
        end
    end
    
end

%% ------------------------函数部分------------------------ %%
% --函数备注：保存数据
function save_xlsx(Path_pro, Name_Tower, name_uncor_cor, name_mean_or_not, Year, Data, Date, Time, SZA, Norm_wl)
    % ---------------------------- 设置输出路径--------------------------- %                                              
    Path_out = strcat(Path_pro, Name_Tower, '_', Year, '_', name_uncor_cor, '_SIF_', name_mean_or_not, '.xlsx');
    % 整理
    wl = Norm_wl;
    if strcmp(name_mean_or_not, 'all')
        Irr_Rad_Refl = [Date, Time, SZA, ...
                    Data(:, 1), Data(:, 2), Data(:, 3), Data(:, 4), Data(:, 5), ...
                    Data(:, 6), Data(:, 7), Data(:, 8), Data(:, 9), Data(:, 10), ...
                    Data(:, 1)./Data(:, 6), Data(:, 2)./Data(:, 7), Data(:, 3)./Data(:, 8), Data(:, 4)./Data(:, 9), Data(:, 5)./Data(:, 10)];
        % 整理SIF数据
        Ha_SIF = [Date, Time, SZA, Data(:, 11), Data(:, 12), Data(:, 13), Data(:, 14), Data(:, 15), Data(:, 16), Data(:, 17)];
        O2B_SIF = [Date, Time, SZA, Data(:, 18), Data(:, 19), Data(:, 20), Data(:, 21), Data(:, 22), Data(:, 23), Data(:, 24)];
        H2O_SIF = [Date, Time, SZA, Data(:, 25), Data(:, 26), Data(:, 27), Data(:, 28), Data(:, 29), Data(:, 30), Data(:, 31)];
        O2A_SIF = [Date, Time, SZA, Data(:, 32), Data(:, 33), Data(:, 34), Data(:, 35), Data(:, 36), Data(:, 37), Data(:, 38)];  
        SIF_fullband = Data(:, 39:end); 

        % ----------------------------------------------输出Excel结果------------------------------------------%
        Names_irr_rad_refl = {'Date', 'Time', 'SZA', 'Rad_680', 'Rad_705', 'Rad_740', 'Rad_770', 'Rad_800', ...
                                             'Irr_680', 'Irr_705', 'Irr_740', 'Irr_770', 'Irr_800', ...
                                             'Refl_680', 'Refl_705', 'Refl_740', 'Refl_770', 'Refl_800'};
        Names_singleband = {'Date', 'Time', 'SZA', 'sFLD', '3FLD', 'iFLD', 'pFLD', 'SFM', 'SVD', 'DOAS'};
        % 输出Irr_Rad_Refl
        xlswrite(Path_out, Names_irr_rad_refl , 'Irr_Rad_Refl', 'A1');
        xlswrite(Path_out, Irr_Rad_Refl , 'Irr_Rad_Refl', 'A2');
        % 输出Ha
        xlswrite(Path_out, Names_singleband , 'SIF_Ha', 'A1');
        xlswrite(Path_out, Ha_SIF , 'SIF_Ha', 'A2');
        % 输出O2B
        xlswrite(Path_out, Names_singleband , 'SIF_O2B', 'A1');
        xlswrite(Path_out, O2B_SIF , 'SIF_O2B', 'A2');
        % 输出H2O
        xlswrite(Path_out, Names_singleband , 'SIF_H2O', 'A1');
        xlswrite(Path_out, H2O_SIF , 'SIF_H2O', 'A2');
        % 输出O2A
        xlswrite(Path_out, Names_singleband , 'SIF_O2A', 'A1');
        xlswrite(Path_out, O2A_SIF , 'SIF_O2A', 'A2');
        % 输出fullband
        xlswrite(Path_out, Date , 'SIF_Fullband', 'A2');
        xlswrite(Path_out, Time , 'SIF_Fullband', 'B2');
        xlswrite(Path_out, wl' , 'SIF_Fullband', 'C1');
        xlswrite(Path_out, SIF_fullband, 'SIF_Fullband', 'C2');
    elseif strcmp(name_mean_or_not, '30min')
        Irr_Rad_Refl = [Date, Time, SZA, ...
                    Data(:, 1), Data(:, 2), Data(:, 3), Data(:, 4), Data(:, 5), ...
                    Data(:, 6), Data(:, 7), Data(:, 8), Data(:, 9), Data(:, 10), ...
                    Data(:, 1)./Data(:, 6), Data(:, 2)./Data(:, 7), Data(:, 3)./Data(:, 8), Data(:, 4)./Data(:, 9), Data(:, 5)./Data(:, 10)];
        % 整理SIF数据
        Ha_SIF = [Date, Time, SZA, Data(:, 11), Data(:, 12), Data(:, 13), Data(:, 14), Data(:, 15), Data(:, 16), Data(:, 17)];
        O2B_SIF = [Date, Time, SZA, Data(:, 18), Data(:, 19), Data(:, 20), Data(:, 21), Data(:, 22), Data(:, 23), Data(:, 24)];
        H2O_SIF = [Date, Time, SZA, Data(:, 25), Data(:, 26), Data(:, 27), Data(:, 28), Data(:, 29), Data(:, 30), Data(:, 31)];
        O2A_SIF = [Date, Time, SZA, Data(:, 32), Data(:, 33), Data(:, 34), Data(:, 35), Data(:, 36), Data(:, 37), Data(:, 38)];  
        SIF_fullband = Data(:, 39:end); 

        % ----------------------------------------------输出Excel结果------------------------------------------%
        Names_irr_rad_refl = {'Date_30min', 'Time_30min', 'SZA_30min', 'Rad_680_30min', 'Rad_705_30min', 'Rad_740_30min', 'Rad_770_30min', 'Rad_800_30min', ...
                                             'Irr_680_30min', 'Irr_705_30min', 'Irr_740_30min', 'Irr_770_30min', 'Irr_800_30min', ...
                                             'Refl_680_30min', 'Refl_705_30min', 'Refl_740_30min', 'Refl_770_30min', 'Refl_800_30min'};
        Names_singleband = {'Date_30min', 'Time_30min', 'SZA_30min', 'sFLD_30min', '3FLD_30min', 'iFLD_30min', 'pFLD_30min', 'SFM_30min', 'SVD_30min', 'DOAS_30min'};
        % 输出Irr_Rad_Refl
        xlswrite(Path_out, Names_irr_rad_refl , 'Irr_Rad_Refl_30min', 'A1');
        xlswrite(Path_out, Irr_Rad_Refl , 'Irr_Rad_Refl_30min', 'A2');
        % 输出Ha
        xlswrite(Path_out, Names_singleband , 'SIF_Ha_30min', 'A1');
        xlswrite(Path_out, Ha_SIF , 'SIF_Ha_30min', 'A2');
        % 输出O2B
        xlswrite(Path_out, Names_singleband , 'SIF_O2B_30min', 'A1');
        xlswrite(Path_out, O2B_SIF , 'SIF_O2B_30min', 'A2');
        % 输出H2O
        xlswrite(Path_out, Names_singleband , 'SIF_H2O_30min', 'A1');
        xlswrite(Path_out, H2O_SIF , 'SIF_H2O_30min', 'A2');
        % 输出O2A
        xlswrite(Path_out, Names_singleband , 'SIF_O2A_30min', 'A1');
        xlswrite(Path_out, O2A_SIF , 'SIF_O2A_30min', 'A2');
        % 输出fullband
        xlswrite(Path_out, Date , 'SIF_Fullband_30min', 'A2');
        xlswrite(Path_out, Time , 'SIF_Fullband_30min', 'B2');
        xlswrite(Path_out, wl' , 'SIF_Fullband_30min', 'C1');
        xlswrite(Path_out, SIF_fullband, 'SIF_Fullband_30min', 'C2');
    end
end
function save_mat(Path_pro, Name_Tower, name_uncor_cor, name_mean_or_not, Year, Data, Date, Time, SZA, Norm_wl)
    % ---------------------------- DEFINE THE FILE --------------------------- %      
    Path_out = strcat(Path_pro, Name_Tower, '_', Year, '_', name_uncor_cor, '_SIF_', name_mean_or_not, '.mat');
    %mW/m2/nm/sr
    wl = Norm_wl;
    if strcmp(name_mean_or_not, 'all')
        Rad_680 = Data(:,1);
        Rad_705 = Data(:,2);
        Rad_740 = Data(:,3);
        Rad_770 = Data(:,4);
        Rad_800 = Data(:,5);
        Irr_680 = Data(:,6);       
        Irr_705 = Data(:,7);
        Irr_740 = Data(:,8);
        Irr_770 = Data(:,9);
        Irr_800 = Data(:,10);
        Refl_680 = Rad_680./Irr_680;
        Refl_705 = Rad_705./Irr_705;
        Refl_740 = Rad_740./Irr_740;
        Refl_770 = Rad_770./Irr_770;
        Refl_800 = Rad_800./Irr_800;

        SIF_Ha_sFLD = Data(:,11);
        SIF_Ha_3FLD= Data(:,12);
        SIF_Ha_iFLD = Data(:,13);
        SIF_Ha_pFLD = Data(:,14);
        SIF_Ha_SFM = Data(:,15);
        SIF_Ha_SVD = Data(:,16);
        SIF_Ha_DOAS = Data(:,17);

        SIF_O2B_sFLD = Data(:,18);
        SIF_O2B_3FLD = Data(:,19);
        SIF_O2B_iFLD = Data(:,20);
        SIF_O2B_pFLD = Data(:,21);
        SIF_O2B_SFM = Data(:,22);
        SIF_O2B_SVD = Data(:,23);
        SIF_O2B_DOAS = Data(:,24);

        SIF_H2O_sFLD = Data(:,25);
        SIF_H2O_3FLD = Data(:,26);
        SIF_H2O_iFLD = Data(:,27);
        SIF_H2O_pFLD = Data(:,28);
        SIF_H2O_SFM = Data(:,29);
        SIF_H2O_SVD = Data(:,30);
        SIF_H2O_DOAS = Data(:,31);

        SIF_O2A_sFLD = Data(:,32);
        SIF_O2A_3FLD = Data(:,33);
        SIF_O2A_iFLD = Data(:,34);
        SIF_O2A_pFLD = Data(:,35);
        SIF_O2A_SFM = Data(:,36);
        SIF_O2A_SVD = Data(:,37);
        SIF_O2A_DOAS = Data(:,38);
        SIF_Fullband = Data(:,39:end);
    
        save(Path_out,...
        'Rad_680', 'Rad_705', 'Rad_740', 'Rad_770', 'Rad_800', ...
        'Irr_680', 'Irr_705', 'Irr_740', 'Irr_770', 'Irr_800', ...
        'Refl_680', 'Refl_705', 'Refl_740', 'Refl_770', 'Refl_800', ... 
        'SIF_Ha_sFLD', 'SIF_Ha_3FLD', 'SIF_Ha_iFLD', 'SIF_Ha_pFLD', 'SIF_Ha_SFM', 'SIF_Ha_SVD', 'SIF_Ha_DOAS', ...
        'SIF_O2B_sFLD', 'SIF_O2B_3FLD', 'SIF_O2B_iFLD', 'SIF_O2B_pFLD', 'SIF_O2B_SFM', 'SIF_O2B_SVD', 'SIF_O2B_DOAS', ...
        'SIF_H2O_sFLD', 'SIF_H2O_3FLD', 'SIF_H2O_iFLD', 'SIF_H2O_pFLD', 'SIF_H2O_SFM', 'SIF_H2O_SVD', 'SIF_H2O_DOAS', ...
        'SIF_O2A_sFLD', 'SIF_O2A_3FLD', 'SIF_O2A_iFLD', 'SIF_O2A_pFLD', 'SIF_O2A_SFM', 'SIF_O2A_SVD', 'SIF_O2A_DOAS', ...
        'SIF_Fullband', ...
        'wl', ...
        'Date', ...
        'Time', ...
        'SZA' );
    
    elseif strcmp(name_mean_or_not, '30min')
        Date_30min = Date;
        Time_30min = Time;
        SZA_30min = SZA;
        Rad_680_30min = Data(:,1);
        Rad_705_30min = Data(:,2);
        Rad_740_30min = Data(:,3);
        Rad_770_30min = Data(:,4);
        Rad_800_30min = Data(:,5);
        Irr_680_30min = Data(:,6);       
        Irr_705_30min = Data(:,7);
        Irr_740_30min = Data(:,8);
        Irr_770_30min = Data(:,9);
        Irr_800_30min = Data(:,10);
        Refl_680_30min=Rad_680_30min./Irr_680_30min;
        Refl_705_30min=Rad_705_30min./Irr_705_30min;
        Refl_740_30min=Rad_740_30min./Irr_740_30min;
        Refl_770_30min=Rad_770_30min./Irr_770_30min;
        Refl_800_30min=Rad_800_30min./Irr_800_30min;

        SIF_Ha_sFLD_30min = Data(:,11);
        SIF_Ha_3FLD_30min = Data(:,12);
        SIF_Ha_iFLD_30min = Data(:,13);
        SIF_Ha_pFLD_30min = Data(:,14);
        SIF_Ha_SFM_30min = Data(:,15);
        SIF_Ha_SVD_30min = Data(:,16);
        SIF_Ha_DOAS_30min = Data(:,17);

        SIF_O2B_sFLD_30min = Data(:,18);
        SIF_O2B_3FLD_30min = Data(:,19);
        SIF_O2B_iFLD_30min = Data(:,20);
        SIF_O2B_pFLD_30min = Data(:,21);
        SIF_O2B_SFM_30min = Data(:,22);
        SIF_O2B_SVD_30min = Data(:,23);
        SIF_O2B_DOAS_30min = Data(:,24);

        SIF_H2O_sFLD_30min = Data(:,25);
        SIF_H2O_3FLD_30min = Data(:,26);
        SIF_H2O_iFLD_30min = Data(:,27);
        SIF_H2O_pFLD_30min = Data(:,28);
        SIF_H2O_SFM_30min = Data(:,29);
        SIF_H2O_SVD_30min = Data(:,30);
        SIF_H2O_DOAS_30min = Data(:,31);

        SIF_O2A_sFLD_30min = Data(:,32);
        SIF_O2A_3FLD_30min = Data(:,33);
        SIF_O2A_iFLD_30min = Data(:,34);
        SIF_O2A_pFLD_30min = Data(:,35);
        SIF_O2A_SFM_30min = Data(:,36);
        SIF_O2A_SVD_30min = Data(:,37);
        SIF_O2A_DOAS_30min = Data(:,38);

        SIF_Fullband_30min = Data(:,39:end);
        
        save(Path_out,...
        'Rad_680_30min', 'Rad_705_30min', 'Rad_740_30min', 'Rad_770_30min', 'Rad_800_30min', ...
        'Irr_680_30min', 'Irr_705_30min', 'Irr_740_30min', 'Irr_770_30min', 'Irr_800_30min', ...
        'Refl_680_30min', 'Refl_705_30min', 'Refl_740_30min', 'Refl_770_30min', 'Refl_800_30min', ... 
        'SIF_Ha_sFLD_30min', 'SIF_Ha_3FLD_30min', 'SIF_Ha_iFLD_30min', 'SIF_Ha_pFLD_30min', 'SIF_Ha_SFM_30min', 'SIF_Ha_SVD_30min', 'SIF_Ha_DOAS_30min', ...
        'SIF_O2B_sFLD_30min', 'SIF_O2B_3FLD_30min', 'SIF_O2B_iFLD_30min', 'SIF_O2B_pFLD_30min', 'SIF_O2B_SFM_30min', 'SIF_O2B_SVD_30min', 'SIF_O2B_DOAS_30min', ...
        'SIF_H2O_sFLD_30min', 'SIF_H2O_3FLD_30min', 'SIF_H2O_iFLD_30min', 'SIF_H2O_pFLD_30min', 'SIF_H2O_SFM_30min', 'SIF_H2O_SVD_30min', 'SIF_H2O_DOAS_30min', ...
        'SIF_O2A_sFLD_30min', 'SIF_O2A_3FLD_30min', 'SIF_O2A_iFLD_30min', 'SIF_O2A_pFLD_30min', 'SIF_O2A_SFM_30min', 'SIF_O2A_SVD_30min', 'SIF_O2A_DOAS_30min', ...
        'SIF_Fullband_30min', ...
        'wl', ...
        'Date_30min', ...
        'Time_30min', ...
        'SZA_30min' );
    end
    
end
function save_nc(Path_pro, Name_Tower, name_uncor_cor, name_mean_or_not, Year, Data, Date, Time, SZA, Norm_wl)
    % ---------------------------- 设置输出路径--------------------------- %                                              
    Path_out = strcat(Path_pro, Name_Tower, '_', Year, '_', name_uncor_cor, '_SIF_', name_mean_or_not, '.nc');
    
    wl = Norm_wl;
    nmeas = size(Time, 1);
    npixels = size(wl, 1);
    if strcmp(name_mean_or_not, 'all')
        ncid=netcdf.create(Path_out,'CLOBBER'); %创建一个存放数据的nc文件
        %-----------------------------define dimension-----------------------------%   
        dimidx = netcdf.defDim(ncid,'Measures', nmeas);    
        dimidy = netcdf.defDim(ncid,'Wavelength',npixels);
        %----------------------------define new variables---------------------------------%
        varid1 = netcdf.defVar(ncid,'Date','double',[dimidx]);
        varid2 = netcdf.defVar(ncid,'Time','double',[dimidx]);
        varid3 = netcdf.defVar(ncid,'SZA','double',[dimidx]);
        varid4 = netcdf.defVar(ncid,'wl','double',[dimidy]);
        varid5 = netcdf.defVar(ncid,'Irr_680','double',[dimidx]);
        varid6 = netcdf.defVar(ncid,'Irr_705','double',[dimidx]);
        varid7 = netcdf.defVar(ncid,'Irr_740','double',[dimidx]);
        varid8 = netcdf.defVar(ncid,'Irr_770','double',[dimidx]);
        varid9 = netcdf.defVar(ncid,'Irr_800','double',[dimidx]);
        varid10 = netcdf.defVar(ncid,'Rad_680','double',[dimidx]);
        varid11 = netcdf.defVar(ncid,'Rad_705','double',[dimidx]);
        varid12 = netcdf.defVar(ncid,'Rad_740','double',[dimidx]);
        varid13 = netcdf.defVar(ncid,'Rad_770','double',[dimidx]);
        varid14 = netcdf.defVar(ncid,'Rad_800','double',[dimidx]);
        varid15 = netcdf.defVar(ncid,'Refl_680','double',[dimidx]);
        varid16 = netcdf.defVar(ncid,'Refl_705','double',[dimidx]);
        varid17 = netcdf.defVar(ncid,'Refl_740','double',[dimidx]);
        varid18 = netcdf.defVar(ncid,'Refl_770','double',[dimidx]);
        varid19 = netcdf.defVar(ncid,'Refl_800','double',[dimidx]);
        % Ha大气吸收窗口
        varid20 = netcdf.defVar(ncid,'SIF_Ha_sFLD','double',[dimidx]);
        varid21 = netcdf.defVar(ncid,'SIF_Ha_3FLD','double',[dimidx]);
        varid22 = netcdf.defVar(ncid,'SIF_Ha_iFLD','double',[dimidx]);
        varid23 = netcdf.defVar(ncid,'SIF_Ha_pFLD','double',[dimidx]);
        varid24 = netcdf.defVar(ncid,'SIF_Ha_SFM','double',[dimidx]);
        varid25 = netcdf.defVar(ncid,'SIF_Ha_SVD','double',[dimidx]);
        varid26 = netcdf.defVar(ncid,'SIF_Ha_DOAS','double',[dimidx]);
        % O2B大气吸收窗口
        varid27 = netcdf.defVar(ncid,'SIF_O2B_sFLD','double',[dimidx]);
        varid28 = netcdf.defVar(ncid,'SIF_O2B_3FLD','double',[dimidx]);
        varid29 = netcdf.defVar(ncid,'SIF_O2B_iFLD','double',[dimidx]);
        varid30 = netcdf.defVar(ncid,'SIF_O2B_pFLD','double',[dimidx]);
        varid31 = netcdf.defVar(ncid,'SIF_O2B_SFM','double',[dimidx]);
        varid32 = netcdf.defVar(ncid,'SIF_O2B_SVD','double',[dimidx]);
        varid33 = netcdf.defVar(ncid,'SIF_O2B_DOAS','double',[dimidx]);
        % H2O大气吸收窗口
        varid34 = netcdf.defVar(ncid,'SIF_H2O_sFLD','double',[dimidx]);
        varid35 = netcdf.defVar(ncid,'SIF_H2O_3FLD','double',[dimidx]);
        varid36 = netcdf.defVar(ncid,'SIF_H2O_iFLD','double',[dimidx]);
        varid37 = netcdf.defVar(ncid,'SIF_H2O_pFLD','double',[dimidx]);
        varid38 = netcdf.defVar(ncid,'SIF_H2O_SFM','double',[dimidx]);
        varid39 = netcdf.defVar(ncid,'SIF_H2O_SVD','double',[dimidx]);
        varid40 = netcdf.defVar(ncid,'SIF_H2O_DOAS','double',[dimidx]);
        % O2A大气吸收窗口
        varid41 = netcdf.defVar(ncid,'SIF_O2A_sFLD','double',[dimidx]);
        varid42 = netcdf.defVar(ncid,'SIF_O2A_3FLD','double',[dimidx]);
        varid43 = netcdf.defVar(ncid,'SIF_O2A_iFLD','double',[dimidx]);
        varid44 = netcdf.defVar(ncid,'SIF_O2A_pFLD','double',[dimidx]);
        varid45 = netcdf.defVar(ncid,'SIF_O2A_SFM','double',[dimidx]);
        varid46 = netcdf.defVar(ncid,'SIF_O2A_SVD','double',[dimidx]);
        varid47 = netcdf.defVar(ncid,'SIF_O2A_DOAS','double',[dimidx]);
        % 全波段反演结果
        varid48 = netcdf.defVar(ncid,'Fullband','double',[dimidx dimidy]);

        %---------------------------define attributes of the new variables--------------%  
        netcdf.putAtt(ncid,varid1,'units','/');                                                     %单位信息和long_name，其它的信息可依此定义
        netcdf.putAtt(ncid,varid2,'units','/');
        netcdf.putAtt(ncid,varid3,'units','degree');      
        netcdf.putAtt(ncid,varid4,'units','nm');  
        netcdf.putAtt(ncid,varid5,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid6,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid7,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid8,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid9,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid10,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid11,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid12,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid13,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid14,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid15,'units','/'); 
        netcdf.putAtt(ncid,varid16,'units','/');
        netcdf.putAtt(ncid,varid17,'units','/');
        netcdf.putAtt(ncid,varid18,'units','/');
        netcdf.putAtt(ncid,varid19,'units','/');
        % Ha大气吸收窗口
        netcdf.putAtt(ncid,varid20,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid21,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid22,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid23,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid24,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid25,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid26,'units','mW/m2/nm/sr');
        % O2B大气吸收窗口
        netcdf.putAtt(ncid,varid27,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid28,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid29,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid30,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid31,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid32,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid33,'units','mW/m2/nm/sr');
        % H2O大气吸收窗口
        netcdf.putAtt(ncid,varid34,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid35,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid36,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid37,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid38,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid39,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid40,'units','mW/m2/nm/sr');
        % O2A大气吸收窗口
        netcdf.putAtt(ncid,varid41,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid42,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid43,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid44,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid45,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid46,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid47,'units','mW/m2/nm/sr');
        % 全波段反演结果
        netcdf.putAtt(ncid,varid48,'units','mW/m2/nm/sr');

        netcdf.putAtt(ncid,varid1,'long_name','Date');                                                     %单位信息和long_name，其它的信息可依此定义
        netcdf.putAtt(ncid,varid2,'long_name','decimal Julian day');
        netcdf.putAtt(ncid,varid3,'long_name','Solar zenith angle');   
        netcdf.putAtt(ncid,varid4,'long_name','Wavelength');  
        netcdf.putAtt(ncid,varid5,'long_name','Solar irradiance / pi at 680 nm [4nm bin]');
        netcdf.putAtt(ncid,varid6,'long_name','Solar irradiance / pi at 705 nm [4nm bin]');
        netcdf.putAtt(ncid,varid7,'long_name','Solar irradiance / pi at 740 nm [4nm bin]');
        netcdf.putAtt(ncid,varid8,'long_name','Solar irradiance / pi at 770 nm [4nm bin]');
        netcdf.putAtt(ncid,varid9,'long_name','Solar irradiance / pi at 800 nm [4nm bin]');
        netcdf.putAtt(ncid,varid10,'long_name','Canopy radiance at 680 nm [4nm bin]');
        netcdf.putAtt(ncid,varid11,'long_name','Canopy radiance at 705 nm [4nm bin]');
        netcdf.putAtt(ncid,varid12,'long_name','Canopy radiance at 740 nm [4nm bin]');
        netcdf.putAtt(ncid,varid13,'long_name','Canopy radiance at 770 nm [4nm bin]');
        netcdf.putAtt(ncid,varid14,'long_name','Canopy radiance at 800 nm [4nm bin]');
        netcdf.putAtt(ncid,varid15,'long_name','Canopy reflectance at 680 nm [4nm bin]');
        netcdf.putAtt(ncid,varid16,'long_name','Canopy reflectance at 705 nm [4nm bin]');
        netcdf.putAtt(ncid,varid17,'long_name','Canopy reflectance at 740 nm [4nm bin]');
        netcdf.putAtt(ncid,varid18,'long_name','Canopy reflectance at 770 nm [4nm bin]');
        netcdf.putAtt(ncid,varid19,'long_name','Canopy reflectance at 800 nm [4nm bin]');
        % Ha大气吸收窗口
        netcdf.putAtt(ncid,varid20,'long_name','Canopy SIF at 657 nm by sFLD');
        netcdf.putAtt(ncid,varid21,'long_name','Canopy SIF at 657 nm by 3FLD');
        netcdf.putAtt(ncid,varid22,'long_name','Canopy SIF at 657 nm by iFLD');
        netcdf.putAtt(ncid,varid23,'long_name','Canopy SIF at 657 nm by pFLD');
        netcdf.putAtt(ncid,varid24,'long_name','Canopy SIF at 657 nm by SFM');
        netcdf.putAtt(ncid,varid25,'long_name','Canopy SIF at 657 nm by SVD');
        netcdf.putAtt(ncid,varid26,'long_name','Canopy SIF at 657 nm by DOAS');
        % O2B大气吸收窗口
        netcdf.putAtt(ncid,varid27,'long_name','Canopy SIF at 687 nm by sFLD');
        netcdf.putAtt(ncid,varid28,'long_name','Canopy SIF at 687 nm by 3FLD');
        netcdf.putAtt(ncid,varid29,'long_name','Canopy SIF at 687 nm by iFLD');
        netcdf.putAtt(ncid,varid30,'long_name','Canopy SIF at 687 nm by pFLD');
        netcdf.putAtt(ncid,varid31,'long_name','Canopy SIF at 687 nm by SFM');
        netcdf.putAtt(ncid,varid32,'long_name','Canopy SIF at 687 nm by SVD');
        netcdf.putAtt(ncid,varid33,'long_name','Canopy SIF at 687 nm by DOAS');
        % H2O大气吸收窗口
        netcdf.putAtt(ncid,varid34,'long_name','Canopy SIF at 719 nm by sFLD');
        netcdf.putAtt(ncid,varid35,'long_name','Canopy SIF at 719 nm by 3FLD');
        netcdf.putAtt(ncid,varid36,'long_name','Canopy SIF at 719 nm by iFLD');
        netcdf.putAtt(ncid,varid37,'long_name','Canopy SIF at 719 nm by pFLD');
        netcdf.putAtt(ncid,varid38,'long_name','Canopy SIF at 719 nm by SFM');
        netcdf.putAtt(ncid,varid39,'long_name','Canopy SIF at 719 nm by SVD');
        netcdf.putAtt(ncid,varid40,'long_name','Canopy SIF at 719 nm by DOAS');
        % O2A大气吸收窗口
        netcdf.putAtt(ncid,varid41,'long_name','Canopy SIF at 761 nm by sFLD');
        netcdf.putAtt(ncid,varid42,'long_name','Canopy SIF at 761 nm by 3FLD');
        netcdf.putAtt(ncid,varid43,'long_name','Canopy SIF at 761 nm by iFLD');
        netcdf.putAtt(ncid,varid44,'long_name','Canopy SIF at 761 nm by pFLD');
        netcdf.putAtt(ncid,varid45,'long_name','Canopy SIF at 761 nm by SFM');
        netcdf.putAtt(ncid,varid46,'long_name','Canopy SIF at 761 nm by SVD');
        netcdf.putAtt(ncid,varid47,'long_name','Canopy SIF at 761 nm by DOAS');
        % 全波段反演结果
        netcdf.putAtt(ncid,varid48,'long_name','Canopy SIF at full-wave band by F_SFM');

        netcdf.endDef(ncid);
        %
        netcdf.putVar(ncid,varid1,Date);
        netcdf.putVar(ncid,varid2,Time);
        netcdf.putVar(ncid,varid3,SZA);
        netcdf.putVar(ncid,varid4,wl);
        netcdf.putVar(ncid,varid5,Data(:, 1));
        netcdf.putVar(ncid,varid6,Data(:, 2));
        netcdf.putVar(ncid,varid7,Data(:, 3));
        netcdf.putVar(ncid,varid8,Data(:, 4));
        netcdf.putVar(ncid,varid9,Data(:, 5));
        netcdf.putVar(ncid,varid10,Data(:, 6));
        netcdf.putVar(ncid,varid11,Data(:, 7));
        netcdf.putVar(ncid,varid12,Data(:, 8));
        netcdf.putVar(ncid,varid13,Data(:, 9));
        netcdf.putVar(ncid,varid14,Data(:, 10));
        netcdf.putVar(ncid,varid15,Data(:, 1)./Data(:, 6));
        netcdf.putVar(ncid,varid16,Data(:, 2)./Data(:, 7));
        netcdf.putVar(ncid,varid17,Data(:, 3)./Data(:, 8));
        netcdf.putVar(ncid,varid18,Data(:, 4)./Data(:, 9));
        netcdf.putVar(ncid,varid19,Data(:, 5)./Data(:, 10));
        % Ha大气吸收窗口
        netcdf.putVar(ncid,varid20,Data(:, 11));
        netcdf.putVar(ncid,varid21,Data(:, 12));
        netcdf.putVar(ncid,varid22,Data(:, 13));
        netcdf.putVar(ncid,varid23,Data(:, 14));
        netcdf.putVar(ncid,varid24,Data(:, 15));
        netcdf.putVar(ncid,varid25,Data(:, 16));
        netcdf.putVar(ncid,varid26,Data(:, 17));
        % O2B大气吸收窗口
        netcdf.putVar(ncid,varid27,Data(:, 18));
        netcdf.putVar(ncid,varid28,Data(:, 19));
        netcdf.putVar(ncid,varid29,Data(:, 20));
        netcdf.putVar(ncid,varid30,Data(:, 21));
        netcdf.putVar(ncid,varid31,Data(:, 22));
        netcdf.putVar(ncid,varid32,Data(:, 23));
        netcdf.putVar(ncid,varid33,Data(:, 24));
        % H2O大气吸收窗口
        netcdf.putVar(ncid,varid34,Data(:, 25));
        netcdf.putVar(ncid,varid35,Data(:, 26));
        netcdf.putVar(ncid,varid36,Data(:, 27));
        netcdf.putVar(ncid,varid37,Data(:, 28));
        netcdf.putVar(ncid,varid38,Data(:, 29));
        netcdf.putVar(ncid,varid39,Data(:, 30));
        netcdf.putVar(ncid,varid40,Data(:, 31));
        % O2A大气吸收窗口
        netcdf.putVar(ncid,varid41,Data(:, 32));
        netcdf.putVar(ncid,varid42,Data(:, 33));
        netcdf.putVar(ncid,varid43,Data(:, 34));
        netcdf.putVar(ncid,varid44,Data(:, 35));
        netcdf.putVar(ncid,varid45,Data(:, 36));
        netcdf.putVar(ncid,varid46,Data(:, 37));
        netcdf.putVar(ncid,varid47,Data(:, 38));
        % 全波段反演结果
        netcdf.putVar(ncid,varid48,Data(:, 39:end));

        netcdf.close(ncid);
    elseif strcmp(name_mean_or_not, '30min')
        ncid=netcdf.create(Path_out,'CLOBBER'); %创建一个存放数据的nc文件
        %-----------------------------define dimension-----------------------------%   
        dimidx = netcdf.defDim(ncid,'Measures_30min', nmeas);    
        dimidy = netcdf.defDim(ncid,'Wavelength',npixels);
        %----------------------------define new variables---------------------------------%
        varid1 = netcdf.defVar(ncid,'Date_30min','double',[dimidx]);
        varid2 = netcdf.defVar(ncid,'Time_30min','double',[dimidx]);
        varid3 = netcdf.defVar(ncid,'SZA_30min','double',[dimidx]);
        varid4 = netcdf.defVar(ncid,'wl','double',[dimidy]);
        varid5 = netcdf.defVar(ncid,'Irr_680_30min','double',[dimidx]);
        varid6 = netcdf.defVar(ncid,'Irr_705_30min','double',[dimidx]);
        varid7 = netcdf.defVar(ncid,'Irr_740_30min','double',[dimidx]);
        varid8 = netcdf.defVar(ncid,'Irr_770_30min','double',[dimidx]);
        varid9 = netcdf.defVar(ncid,'Irr_800_30min','double',[dimidx]);
        varid10 = netcdf.defVar(ncid,'Rad_680_30min','double',[dimidx]);
        varid11 = netcdf.defVar(ncid,'Rad_705_30min','double',[dimidx]);
        varid12 = netcdf.defVar(ncid,'Rad_740_30min','double',[dimidx]);
        varid13 = netcdf.defVar(ncid,'Rad_770_30min','double',[dimidx]);
        varid14 = netcdf.defVar(ncid,'Rad_800_30min','double',[dimidx]);
        varid15 = netcdf.defVar(ncid,'Refl_680_30min','double',[dimidx]);
        varid16 = netcdf.defVar(ncid,'Refl_705_30min','double',[dimidx]);
        varid17 = netcdf.defVar(ncid,'Refl_740_30min','double',[dimidx]);
        varid18 = netcdf.defVar(ncid,'Refl_770_30min','double',[dimidx]);
        varid19 = netcdf.defVar(ncid,'Refl_800_30min','double',[dimidx]);
        % Ha大气吸收窗口
        varid20 = netcdf.defVar(ncid,'SIF_Ha_sFLD_30min','double',[dimidx]);
        varid21 = netcdf.defVar(ncid,'SIF_Ha_3FLD_30min','double',[dimidx]);
        varid22 = netcdf.defVar(ncid,'SIF_Ha_iFLD_30min','double',[dimidx]);
        varid23 = netcdf.defVar(ncid,'SIF_Ha_pFLD_30min','double',[dimidx]);
        varid24 = netcdf.defVar(ncid,'SIF_Ha_SFM_30min','double',[dimidx]);
        varid25 = netcdf.defVar(ncid,'SIF_Ha_SVD_30min','double',[dimidx]);
        varid26 = netcdf.defVar(ncid,'SIF_Ha_DOAS_30min','double',[dimidx]);
        % O2B大气吸收窗口
        varid27 = netcdf.defVar(ncid,'SIF_O2B_sFLD_30min','double',[dimidx]);
        varid28 = netcdf.defVar(ncid,'SIF_O2B_3FLD_30min','double',[dimidx]);
        varid29 = netcdf.defVar(ncid,'SIF_O2B_iFLD_30min','double',[dimidx]);
        varid30 = netcdf.defVar(ncid,'SIF_O2B_pFLD_30min','double',[dimidx]);
        varid31 = netcdf.defVar(ncid,'SIF_O2B_SFM_30min','double',[dimidx]);
        varid32 = netcdf.defVar(ncid,'SIF_O2B_SVD_30min','double',[dimidx]);
        varid33 = netcdf.defVar(ncid,'SIF_O2B_DOAS_30min','double',[dimidx]);
        % H2O大气吸收窗口
        varid34 = netcdf.defVar(ncid,'SIF_H2O_sFLD_30min','double',[dimidx]);
        varid35 = netcdf.defVar(ncid,'SIF_H2O_3FLD_30min','double',[dimidx]);
        varid36 = netcdf.defVar(ncid,'SIF_H2O_iFLD_30min','double',[dimidx]);
        varid37 = netcdf.defVar(ncid,'SIF_H2O_pFLD_30min','double',[dimidx]);
        varid38 = netcdf.defVar(ncid,'SIF_H2O_SFM_30min','double',[dimidx]);
        varid39 = netcdf.defVar(ncid,'SIF_H2O_SVD_30min','double',[dimidx]);
        varid40 = netcdf.defVar(ncid,'SIF_H2O_DOAS_30min','double',[dimidx]);
        % O2A大气吸收窗口
        varid41 = netcdf.defVar(ncid,'SIF_O2A_sFLD_30min','double',[dimidx]);
        varid42 = netcdf.defVar(ncid,'SIF_O2A_3FLD_30min','double',[dimidx]);
        varid43 = netcdf.defVar(ncid,'SIF_O2A_iFLD_30min','double',[dimidx]);
        varid44 = netcdf.defVar(ncid,'SIF_O2A_pFLD_30min','double',[dimidx]);
        varid45 = netcdf.defVar(ncid,'SIF_O2A_SFM_30min','double',[dimidx]);
        varid46 = netcdf.defVar(ncid,'SIF_O2A_SVD_30min','double',[dimidx]);
        varid47 = netcdf.defVar(ncid,'SIF_O2A_DOAS_30min','double',[dimidx]);
        % 全波段反演结果
        varid48 = netcdf.defVar(ncid,'SIF_Fullband_30min','double',[dimidx dimidy]);

        %---------------------------define attributes of the new variables--------------%  
        netcdf.putAtt(ncid,varid1,'units','/');                                                     %单位信息和long_name，其它的信息可依此定义
        netcdf.putAtt(ncid,varid2,'units','/');
        netcdf.putAtt(ncid,varid3,'units','degree');      
        netcdf.putAtt(ncid,varid4,'units','nm');  
        netcdf.putAtt(ncid,varid5,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid6,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid7,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid8,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid9,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid10,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid11,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid12,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid13,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid14,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid15,'units','/'); 
        netcdf.putAtt(ncid,varid16,'units','/');
        netcdf.putAtt(ncid,varid17,'units','/');
        netcdf.putAtt(ncid,varid18,'units','/');
        netcdf.putAtt(ncid,varid19,'units','/');
        % Ha大气吸收窗口
        netcdf.putAtt(ncid,varid20,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid21,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid22,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid23,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid24,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid25,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid26,'units','mW/m2/nm/sr');
        % O2B大气吸收窗口
        netcdf.putAtt(ncid,varid27,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid28,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid29,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid30,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid31,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid32,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid33,'units','mW/m2/nm/sr');
        % H2O大气吸收窗口
        netcdf.putAtt(ncid,varid34,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid35,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid36,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid37,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid38,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid39,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid40,'units','mW/m2/nm/sr');
        % O2A大气吸收窗口
        netcdf.putAtt(ncid,varid41,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid42,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid43,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid44,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid45,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid46,'units','mW/m2/nm/sr');
        netcdf.putAtt(ncid,varid47,'units','mW/m2/nm/sr');
        % 全波段反演结果
        netcdf.putAtt(ncid,varid48,'units','mW/m2/nm/sr');

        netcdf.putAtt(ncid,varid1,'long_name','Date');                                                     %单位信息和long_name，其它的信息可依此定义
        netcdf.putAtt(ncid,varid2,'long_name','decimal Julian day');
        netcdf.putAtt(ncid,varid3,'long_name','Solar zenith angle');   
        netcdf.putAtt(ncid,varid4,'long_name','Wavelength');  
        netcdf.putAtt(ncid,varid5,'long_name','Solar irradiance / pi at 680 nm [4nm bin]');
        netcdf.putAtt(ncid,varid6,'long_name','Solar irradiance / pi at 705 nm [4nm bin]');
        netcdf.putAtt(ncid,varid7,'long_name','Solar irradiance / pi at 740 nm [4nm bin]');
        netcdf.putAtt(ncid,varid8,'long_name','Solar irradiance / pi at 770 nm [4nm bin]');
        netcdf.putAtt(ncid,varid9,'long_name','Solar irradiance / pi at 800 nm [4nm bin]');
        netcdf.putAtt(ncid,varid10,'long_name','Canopy radiance at 680 nm [4nm bin]');
        netcdf.putAtt(ncid,varid11,'long_name','Canopy radiance at 705 nm [4nm bin]');
        netcdf.putAtt(ncid,varid12,'long_name','Canopy radiance at 740 nm [4nm bin]');
        netcdf.putAtt(ncid,varid13,'long_name','Canopy radiance at 770 nm [4nm bin]');
        netcdf.putAtt(ncid,varid14,'long_name','Canopy radiance at 800 nm [4nm bin]');
        netcdf.putAtt(ncid,varid15,'long_name','Canopy reflectance at 680 nm [4nm bin]');
        netcdf.putAtt(ncid,varid16,'long_name','Canopy reflectance at 705 nm [4nm bin]');
        netcdf.putAtt(ncid,varid17,'long_name','Canopy reflectance at 740 nm [4nm bin]');
        netcdf.putAtt(ncid,varid18,'long_name','Canopy reflectance at 770 nm [4nm bin]');
        netcdf.putAtt(ncid,varid19,'long_name','Canopy reflectance at 800 nm [4nm bin]');
        % Ha大气吸收窗口
        netcdf.putAtt(ncid,varid20,'long_name','Canopy SIF at 657 nm by sFLD');
        netcdf.putAtt(ncid,varid21,'long_name','Canopy SIF at 657 nm by 3FLD');
        netcdf.putAtt(ncid,varid22,'long_name','Canopy SIF at 657 nm by iFLD');
        netcdf.putAtt(ncid,varid23,'long_name','Canopy SIF at 657 nm by pFLD');
        netcdf.putAtt(ncid,varid24,'long_name','Canopy SIF at 657 nm by SFM');
        netcdf.putAtt(ncid,varid25,'long_name','Canopy SIF at 657 nm by SVD');
        netcdf.putAtt(ncid,varid26,'long_name','Canopy SIF at 657 nm by DOAS');
        % O2B大气吸收窗口
        netcdf.putAtt(ncid,varid27,'long_name','Canopy SIF at 687 nm by sFLD');
        netcdf.putAtt(ncid,varid28,'long_name','Canopy SIF at 687 nm by 3FLD');
        netcdf.putAtt(ncid,varid29,'long_name','Canopy SIF at 687 nm by iFLD');
        netcdf.putAtt(ncid,varid30,'long_name','Canopy SIF at 687 nm by pFLD');
        netcdf.putAtt(ncid,varid31,'long_name','Canopy SIF at 687 nm by SFM');
        netcdf.putAtt(ncid,varid32,'long_name','Canopy SIF at 687 nm by SVD');
        netcdf.putAtt(ncid,varid33,'long_name','Canopy SIF at 687 nm by DOAS');
        % H2O大气吸收窗口
        netcdf.putAtt(ncid,varid34,'long_name','Canopy SIF at 719 nm by sFLD');
        netcdf.putAtt(ncid,varid35,'long_name','Canopy SIF at 719 nm by 3FLD');
        netcdf.putAtt(ncid,varid36,'long_name','Canopy SIF at 719 nm by iFLD');
        netcdf.putAtt(ncid,varid37,'long_name','Canopy SIF at 719 nm by pFLD');
        netcdf.putAtt(ncid,varid38,'long_name','Canopy SIF at 719 nm by SFM');
        netcdf.putAtt(ncid,varid39,'long_name','Canopy SIF at 719 nm by SVD');
        netcdf.putAtt(ncid,varid40,'long_name','Canopy SIF at 719 nm by DOAS');
        % O2A大气吸收窗口
        netcdf.putAtt(ncid,varid41,'long_name','Canopy SIF at 761 nm by sFLD');
        netcdf.putAtt(ncid,varid42,'long_name','Canopy SIF at 761 nm by 3FLD');
        netcdf.putAtt(ncid,varid43,'long_name','Canopy SIF at 761 nm by iFLD');
        netcdf.putAtt(ncid,varid44,'long_name','Canopy SIF at 761 nm by pFLD');
        netcdf.putAtt(ncid,varid45,'long_name','Canopy SIF at 761 nm by SFM');
        netcdf.putAtt(ncid,varid46,'long_name','Canopy SIF at 761 nm by SVD');
        netcdf.putAtt(ncid,varid47,'long_name','Canopy SIF at 761 nm by DOAS');
        % 全波段反演结果
        netcdf.putAtt(ncid,varid48,'long_name','Canopy SIF at full-wave band by F_SFM');

        netcdf.endDef(ncid);
        %
        netcdf.putVar(ncid,varid1,Date);
        netcdf.putVar(ncid,varid2,Time);
        netcdf.putVar(ncid,varid3,SZA);
        netcdf.putVar(ncid,varid4,wl);
        netcdf.putVar(ncid,varid5,Data(:, 1));
        netcdf.putVar(ncid,varid6,Data(:, 2));
        netcdf.putVar(ncid,varid7,Data(:, 3));
        netcdf.putVar(ncid,varid8,Data(:, 4));
        netcdf.putVar(ncid,varid9,Data(:, 5));
        netcdf.putVar(ncid,varid10,Data(:, 6));
        netcdf.putVar(ncid,varid11,Data(:, 7));
        netcdf.putVar(ncid,varid12,Data(:, 8));
        netcdf.putVar(ncid,varid13,Data(:, 9));
        netcdf.putVar(ncid,varid14,Data(:, 10));
        netcdf.putVar(ncid,varid15,Data(:, 1)./Data(:, 6));
        netcdf.putVar(ncid,varid16,Data(:, 2)./Data(:, 7));
        netcdf.putVar(ncid,varid17,Data(:, 3)./Data(:, 8));
        netcdf.putVar(ncid,varid18,Data(:, 4)./Data(:, 9));
        netcdf.putVar(ncid,varid19,Data(:, 5)./Data(:, 10));
        % Ha大气吸收窗口
        netcdf.putVar(ncid,varid20,Data(:, 11));
        netcdf.putVar(ncid,varid21,Data(:, 12));
        netcdf.putVar(ncid,varid22,Data(:, 13));
        netcdf.putVar(ncid,varid23,Data(:, 14));
        netcdf.putVar(ncid,varid24,Data(:, 15));
        netcdf.putVar(ncid,varid25,Data(:, 16));
        netcdf.putVar(ncid,varid26,Data(:, 17));
        % O2B大气吸收窗口
        netcdf.putVar(ncid,varid27,Data(:, 18));
        netcdf.putVar(ncid,varid28,Data(:, 19));
        netcdf.putVar(ncid,varid29,Data(:, 20));
        netcdf.putVar(ncid,varid30,Data(:, 21));
        netcdf.putVar(ncid,varid31,Data(:, 22));
        netcdf.putVar(ncid,varid32,Data(:, 23));
        netcdf.putVar(ncid,varid33,Data(:, 24));
        % H2O大气吸收窗口
        netcdf.putVar(ncid,varid34,Data(:, 25));
        netcdf.putVar(ncid,varid35,Data(:, 26));
        netcdf.putVar(ncid,varid36,Data(:, 27));
        netcdf.putVar(ncid,varid37,Data(:, 28));
        netcdf.putVar(ncid,varid38,Data(:, 29));
        netcdf.putVar(ncid,varid39,Data(:, 30));
        netcdf.putVar(ncid,varid40,Data(:, 31));
        % O2A大气吸收窗口
        netcdf.putVar(ncid,varid41,Data(:, 32));
        netcdf.putVar(ncid,varid42,Data(:, 33));
        netcdf.putVar(ncid,varid43,Data(:, 34));
        netcdf.putVar(ncid,varid44,Data(:, 35));
        netcdf.putVar(ncid,varid45,Data(:, 36));
        netcdf.putVar(ncid,varid46,Data(:, 37));
        netcdf.putVar(ncid,varid47,Data(:, 38));
        % 全波段反演结果
        netcdf.putVar(ncid,varid48,Data(:, 39:end));

        netcdf.close(ncid);
    end
end
