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
% --载入反演窗口
    [windows] = Setup_Window();
% --路径设置
    Path_wl = '.\Settings\Simulated_Datasets\test_2015-07-18-1311_train_2880\wl.dat';               %%修改：
    Path_Ref = '.\Settings\Simulated_Datasets\test_2015-07-18-1311_train_2880\reflectance.dat';     %%修改：
    Path_Fluor = '.\Settings\Simulated_Datasets\test_2015-07-18-1311_train_2880\fluorescence.dat';  %%修改：
    Norm_wl(:, 1) = importdata('.\Settings\wl_pro.txt');  % 修改：标准波长路径
    Scope = [];                                           % 加载模拟数据，并得出主成分，放到Scope中
    [PC8_Ref, PC5_Fluor] = load_scope(Path_wl, Path_Ref, Path_Fluor, Norm_wl);
    Scope.PC8_Ref = PC8_Ref;
    Scope.PC5_Fluor = PC5_Fluor;
    % 可修改：荧光模拟形状数据路径，for DOAS and SVD
    hf = load('.\Settings\Prescribed_shape\hf.mat');    % 根据观测数据的波长对hf进行重采样
    
    uncor_cor = char('Uncor', 'Cor');      % 对大气校正前后的辐亮度数据，进行SIF计算；
    % 可修改：未校正辐亮度数据路径
    Path_Uncor_Cor = strcat('.\Datasets\', Name_Tower, '\处理结果\', Year, '\', Name_Tower, '_Uncor_Cor_Rad_Ref\');
    % 可修改：结果存放路径,未校正和校正
    Path_pro = strcat('.\Datasets\', Name_Tower, '\处理结果\', Year, '\', Name_Tower, '_Uncor_Cor_SIF_Irr_Rad_Refl\');   
    
 % --设置哪种反演方法，0 = off; 1 = on;
    models = [];
    models.sFLD = 1;    
    models.FLD3 = 1;
    models.iFLD = 1;
    models.pFLD = 1;
    models.SFM  = 1;
    models.F_SFM= 1;
    models.SVD  = 1;
    models.DOAS = 1;
 
%% ------读取未校正辐亮度数据、以及气压温度和透过率数据，之后反演SIF------ %%
for i = 1:1:2
    name_uncor_cor = uncor_cor(i, :);
    Path_raw = Path_Uncor_Cor;
    folder = Path_pro;           
    if exist(folder, 'dir')==0    
        mkdir(folder);            
    end
    Path_pro = folder; 
    
    % ----获取未校正和校正后的辐亮度数据，并大反演SIF---- %
    files = dir(strcat(Path_raw, Name_Tower, '_', name_uncor_cor, '*.xlsx'));   % 这里只读取了xlsx文件
    num_files = size(files, 1);                              % xlsx文件数量
    for j = 1:1:num_files
        Path_j = strcat(Path_raw, files(j).name);    % 文件路径
        [~, name, ext] = fileparts(files(j).name);    % 获取指定文件的路径、文件名和扩展名（文件格式）。
        Date = name(isstrprop(name, 'digit'));        % 获取字符串中的数字，即观测日期；strrep也可以实现  
        year_j = str2double(Date(1:4));
        month_j = str2double(Date(5:6));
        day_j = str2double(Date(7:8));
        DOY = day(datetime(year_j, month_j, day_j), 'dayofyear');
        
        disp('Starting： 读取数据' + string(Date) + ', 反演SIF，之后保存');   tic;
        data = [];
        Data_xlsx = xlsread(Path_j, 2);      
        num_meas = (size(Data_xlsx, 2)-1)/2;
        num_wl = size(Data_xlsx, 1)-2;
        wl = Data_xlsx(3:num_wl+2,1);
        time = Data_xlsx(1, 2:num_meas+1);
        sza = Data_xlsx(2, 2:num_meas+1);
        veg = Data_xlsx(3:num_wl+2, 2:num_meas+1);
        sky = Data_xlsx(3:num_wl+2, num_meas+2:2*num_meas+1);
        % 重新算的SZA，也可以用QE65_sza, 注意时区Timezone
        % SZA = Cal_SZA(QE65_time, year_j, month_j, day_j, lon, lat, 8);   
        data.fracDOY = (DOY + time)';       % decimal Julian day, contained in dataset struct
        data.SZA = sza';                    % Zenith angle 
        data.wl = wl;                      % wavelength associated with each pixel of QE Pro. this should exclude the first and last 10 pixels which are optically inactive
        data.sky = sky;                    % incoming irradiance (E), contained in dataset struct. this is size (nmeas, pixels) where pixels excludes optically inactive pixels
        data.veg = veg;                    % outgoing irradiance (L), see above
        data.hf = hf.hf;                   % prescribed shape of fluorescence emission for DOAS and SVD. this is size (pixels,1)
        data.nmeas = size(data.veg, 2);    % number of measurements
        data.npixels = size(data.wl, 1);   % number of pixels

        %---反演SIF,得到结果
        [result] = SIFretrieval(data, models, Scope, windows);
        %---计算red，rededge，near_red的平均值
        index_680 = data.wl > windows.wl_680(1) & data.wl < windows.wl_680(2);
        index_705 = data.wl > windows.wl_705(1) & data.wl < windows.wl_705(2);
        index_740 = data.wl > windows.wl_740(1) & data.wl < windows.wl_740(2);
        index_770 = data.wl > windows.wl_770(1) & data.wl < windows.wl_770(2);
        index_800 = data.wl > windows.wl_800(1) & data.wl < windows.wl_800(2);

        data.rad_sky_680=mean(data.sky(index_680, :),1,'omitnan')';
        data.rad_sky_705=mean(data.sky(index_705, :),1,'omitnan')';
        data.rad_sky_740=mean(data.sky(index_740, :),1,'omitnan')';
        data.rad_sky_770=mean(data.sky(index_770, :),1,'omitnan')';
        data.rad_sky_800=mean(data.sky(index_800, :),1,'omitnan')';

        data.rad_veg_680=mean(data.veg(index_680, :),1,'omitnan')';
        data.rad_veg_705=mean(data.veg(index_705, :),1,'omitnan')';
        data.rad_veg_740=mean(data.veg(index_740, :),1,'omitnan')';
        data.rad_veg_770=mean(data.veg(index_770, :),1,'omitnan')';
        data.rad_veg_800=mean(data.veg(index_800, :),1,'omitnan')';

        data.refl_680=data.rad_veg_680 ./ data.rad_sky_680;
        data.refl_705=data.rad_veg_705 ./ data.rad_sky_705;
        data.refl_740=data.rad_veg_740 ./ data.rad_sky_740;
        data.refl_770=data.rad_veg_770 ./ data.rad_sky_770;
        data.refl_800=data.rad_veg_800 ./ data.rad_sky_800;
        
        
        % --保存数据为nc, xlsx和mat格式
        save_xlsx(Path_pro, Name_Tower, name_uncor_cor, Date, data, result);
        save_mat(Path_pro, Name_Tower, name_uncor_cor, Date, data, result);
        save_nc(Path_pro, Name_Tower, name_uncor_cor, Date, data, result);

    end
end
disp('Time delays: ' + string(toc) + 's');

%% ------------------------函数部分------------------------ %%
% --函数备注：加载Scope模拟数据
function [PC8_Ref, PC5_Fluor] = load_scope(Path_wl, Path_Ref, Path_Fluor, Norm_wl)
    Wl_SCOPE = dlmread(Path_wl, '', [2,240,2,450]);
    Wl_SCOPE = Wl_SCOPE';
    Ref_SCOPE = dlmread(Path_Ref, '', [2,240,2881,450]);
    Ref_SCOPE = Ref_SCOPE';
    Fluor_SCOPE = dlmread(Path_Fluor, '', 2, 0);
    Fluor_SCOPE = Fluor_SCOPE';
    Ref_SCOPE = interp1(Wl_SCOPE(:,1), Ref_SCOPE(:,:), Norm_wl, 'spline');     % ---首先标准波段插值，scope模拟数据--- %
    Fluor_SCOPE = interp1(Wl_SCOPE(:,1), Fluor_SCOPE(:,:), Norm_wl, 'spline');
    % ---这里注意：直射光和散射光O2-A波段吸收深度不同, 插值消除这种不平滑--- %
    WindowA = [758, 770];
    Flag_A  = Norm_wl > WindowA(1) & Norm_wl < WindowA(2);
    Wl0 = Norm_wl;   Wl0(Flag_A, :) = [];   
    Ref0 = Ref_SCOPE;
    Ref0(Flag_A, :) = [];
    Ref_sp = interp1(Wl0, Ref0, Norm_wl, 'spline');
    Ref_SCOPE = Ref_sp;
    % ---获取主成分--- %
    [PCs_Ref, ~, ~] = Cal_PCA(Ref_SCOPE);
    [PCs_Fluor, ~, ~] = Cal_PCA(Fluor_SCOPE);
    PC8_Ref = PCs_Ref(:,1:7);
    PC8_Ref(:,8) = 0;
    PC5_Fluor = PCs_Fluor(:,1:5);
end

% ---函数备注：保存数据
function save_xlsx(Path_pro, Name_Tower, name_uncor_cor, Date, data, result)
    % ---------------------------- 设置输出路径--------------------------- %                                              
    Path_output = strcat(Path_pro, Name_Tower, '_', name_uncor_cor, '_SIF_Irr_Rad_Refl_',Date,'.xlsx');
    % 整理辐亮度数据
    Irr_Rad_Refl = [data.fracDOY, data.SZA, ...
                data.rad_veg_680, data.rad_veg_705, data.rad_veg_740, data.rad_veg_770, data.rad_veg_800, ...
                data.rad_sky_680, data.rad_sky_705, data.rad_sky_740, data.rad_sky_770, data.rad_sky_800, ...
                data.refl_680, data.refl_705, data.refl_740, data.refl_770, data.refl_800];
    % 整理SIF数据
    Ha_SIF = [data.fracDOY, data.SZA, result.Ha.sFLD_SIF, result.Ha.FLD3_SIF, result.Ha.iFLD_SIF, result.Ha.pFLD_SIF, result.Ha.SFM_SIF, result.Ha.SVD_SIF, result.Ha.DOAS_SIF];
    O2B_SIF = [data.fracDOY, data.SZA, result.O2B.sFLD_SIF, result.O2B.FLD3_SIF, result.O2B.iFLD_SIF, result.O2B.pFLD_SIF, result.O2B.SFM_SIF, result.O2B.SVD_SIF, result.O2B.DOAS_SIF];
    H2O_SIF = [data.fracDOY, data.SZA, result.H2O.sFLD_SIF, result.H2O.FLD3_SIF, result.H2O.iFLD_SIF, result.H2O.pFLD_SIF, result.H2O.SFM_SIF, result.H2O.SVD_SIF, result.H2O.DOAS_SIF];
    O2A_SIF = [data.fracDOY, data.SZA, result.O2A.sFLD_SIF, result.O2A.FLD3_SIF, result.O2A.iFLD_SIF, result.O2A.pFLD_SIF, result.O2A.SFM_SIF, result.O2A.SVD_SIF, result.O2A.DOAS_SIF];  
    SIF_fullband = result.FSFM.FSFM_SIF; 
    
    % ----------------------------------------------输出Excel结果------------------------------------------%
    Names_irr_rad_refl = {'Time', 'SZA', 'Rad_680', 'Rad_705', 'Rad_740', 'Rad_770', 'Rad_800', ...
                                         'Irr_680', 'Irr_705', 'Irr_740', 'Irr_770', 'Irr_800', ...
                                         'Refl_680', 'Refl_705', 'Refl_740', 'Refl_770', 'Refl_800'};
    Names_singleband = {'Time', 'SZA', 'sFLD', '3FLD', 'iFLD', 'pFLD', 'SFM', 'SVD', 'DOAS'};
    % 输出Irr_Rad_Refl
    xlswrite(Path_output, Names_irr_rad_refl , 'Irr_Rad_Refl', 'A1');
    xlswrite(Path_output, Irr_Rad_Refl , 'Irr_Rad_Refl', 'A2');
    % 输出Ha
    xlswrite(Path_output, Names_singleband , 'SIF_Ha', 'A1');
    xlswrite(Path_output, Ha_SIF , 'SIF_Ha', 'A2');
    % 输出O2B
    xlswrite(Path_output, Names_singleband , 'SIF_O2B', 'A1');
    xlswrite(Path_output, O2B_SIF , 'SIF_O2B', 'A2');
    % 输出H2O
    xlswrite(Path_output, Names_singleband , 'SIF_H2O', 'A1');
    xlswrite(Path_output, H2O_SIF , 'SIF_H2O', 'A2');
    % 输出O2A
    xlswrite(Path_output, Names_singleband , 'SIF_O2A', 'A1');
    xlswrite(Path_output, O2A_SIF , 'SIF_O2A', 'A2');
    % 输出fullband
    xlswrite(Path_output, data.fracDOY , 'SIF_Fullband', 'A2');
    xlswrite(Path_output, data.wl' , 'SIF_Fullband', 'B1');
    xlswrite(Path_output, SIF_fullband, 'SIF_Fullband', 'B2');
end
function save_mat(Path_pro, Name_Tower, name_uncor_cor, Date, data, result)
    % ---------------------------- 设置输出路径--------------------------- %                                              
    Path_output = strcat(Path_pro, Name_Tower, '_', name_uncor_cor, '_SIF_Irr_Rad_Refl_',Date,'.mat');
    % 时间，sza, wavelength
    Time = data.fracDOY;
    SZA = data.SZA;
    wl = data.wl;
    % 整理辐亮度数据
    Rad_680 = data.rad_veg_680; Rad_705 = data.rad_veg_705; Rad_740 = data.rad_veg_740; Rad_770 = data.rad_veg_770; Rad_800 = data.rad_veg_800;
    Irr_680 = data.rad_sky_680; Irr_705 = data.rad_sky_705; Irr_740 = data.rad_sky_740; Irr_770 = data.rad_sky_770; Irr_800 = data.rad_sky_800;
    Refl_680 = data.refl_680;   Refl_705 = data.refl_705;   Refl_740 = data.refl_740;   Refl_770 = data.refl_770;   Refl_800 = data.refl_800;
    
    % 整理Ha
    SIF_Ha_sFLD = result.Ha.sFLD_SIF;  
    SIF_Ha_3FLD = result.Ha.FLD3_SIF;
    SIF_Ha_iFLD = result.Ha.iFLD_SIF;
    SIF_Ha_pFLD = result.Ha.pFLD_SIF;
    SIF_Ha_SFM = result.Ha.SFM_SIF;
    SIF_Ha_SVD = result.Ha.SVD_SIF;
    SIF_Ha_DOAS = result.Ha.DOAS_SIF;
    % 整理O2B
    SIF_O2B_sFLD = result.O2B.sFLD_SIF;  
    SIF_O2B_3FLD = result.O2B.FLD3_SIF;
    SIF_O2B_iFLD = result.O2B.iFLD_SIF;
    SIF_O2B_pFLD = result.O2B.pFLD_SIF;
    SIF_O2B_SFM = result.O2B.SFM_SIF;
    SIF_O2B_SVD = result.O2B.SVD_SIF;
    SIF_O2B_DOAS = result.O2B.DOAS_SIF;
    % 整理H2O
    SIF_H2O_sFLD = result.H2O.sFLD_SIF;  
    SIF_H2O_3FLD = result.H2O.FLD3_SIF;
    SIF_H2O_iFLD = result.H2O.iFLD_SIF;
    SIF_H2O_pFLD = result.H2O.pFLD_SIF;
    SIF_H2O_SFM = result.H2O.SFM_SIF;
    SIF_H2O_SVD = result.H2O.SVD_SIF;
    SIF_H2O_DOAS = result.H2O.DOAS_SIF;
    % 整理O2A
    SIF_O2A_sFLD = result.O2A.sFLD_SIF;  
    SIF_O2A_3FLD = result.O2A.FLD3_SIF;
    SIF_O2A_iFLD = result.O2A.iFLD_SIF;
    SIF_O2A_pFLD = result.O2A.pFLD_SIF;
    SIF_O2A_SFM = result.O2A.SFM_SIF;
    SIF_O2A_SVD = result.O2A.SVD_SIF;
    SIF_O2A_DOAS = result.O2A.DOAS_SIF;
    % 整理全波段
    SIF_Fullband = result.FSFM.FSFM_SIF;
    
    % ---------------------输出mat结果-----------------------%
    save(Path_output, 'Time', 'SZA', 'wl', ...
                      'Rad_680', 'Rad_705', 'Rad_740', 'Rad_770', 'Rad_800', ...
                      'Irr_680', 'Irr_705', 'Irr_740', 'Irr_770', 'Irr_800', ...
                      'Refl_680', 'Refl_705', 'Refl_740', 'Refl_770', 'Refl_800', ...
         'SIF_Ha_sFLD', 'SIF_Ha_3FLD', 'SIF_Ha_iFLD', 'SIF_Ha_pFLD', 'SIF_Ha_SFM', 'SIF_Ha_SVD', 'SIF_Ha_DOAS', ...
         'SIF_O2B_sFLD', 'SIF_O2B_3FLD', 'SIF_O2B_iFLD', 'SIF_O2B_pFLD', 'SIF_O2B_SFM', 'SIF_O2B_SVD', 'SIF_O2B_DOAS', ...
         'SIF_H2O_sFLD', 'SIF_H2O_3FLD', 'SIF_H2O_iFLD', 'SIF_H2O_pFLD', 'SIF_H2O_SFM', 'SIF_H2O_SVD', 'SIF_H2O_DOAS', ...
         'SIF_O2A_sFLD', 'SIF_O2A_3FLD', 'SIF_O2A_iFLD', 'SIF_O2A_pFLD', 'SIF_O2A_SFM', 'SIF_O2A_SVD', 'SIF_O2A_DOAS', ...
         'SIF_Fullband');
end
function save_nc(Path_pro, Name_Tower, name_uncor_cor, Date, data, result)
    % ---------------------------- 设置输出路径--------------------------- %                                              
    Path_output = strcat(Path_pro, Name_Tower, '_', name_uncor_cor, '_SIF_Irr_Rad_Refl_',Date,'.nc');
    ncid=netcdf.create(Path_output,'CLOBBER'); %创建一个存放数据的nc文件
    %-----------------------------define dimension-----------------------------%   
    dimidx = netcdf.defDim(ncid,'Measures', data.nmeas);    
    dimidy = netcdf.defDim(ncid,'Wavelength',data.npixels);
    %----------------------------define new variables---------------------------------%
    varid1 = netcdf.defVar(ncid,'Time','double',[dimidx]);
    varid2 = netcdf.defVar(ncid,'SZA','double',[dimidx]);
    varid3 = netcdf.defVar(ncid,'wl','double',[dimidy]);
    varid4 = netcdf.defVar(ncid,'Irr_680','double',[dimidx]);
    varid5 = netcdf.defVar(ncid,'Irr_705','double',[dimidx]);
    varid6 = netcdf.defVar(ncid,'Irr_740','double',[dimidx]);
    varid7 = netcdf.defVar(ncid,'Irr_770','double',[dimidx]);
    varid8 = netcdf.defVar(ncid,'Irr_800','double',[dimidx]);
    varid9 = netcdf.defVar(ncid,'Rad_680','double',[dimidx]);
    varid10 = netcdf.defVar(ncid,'Rad_705','double',[dimidx]);
    varid11 = netcdf.defVar(ncid,'Rad_740','double',[dimidx]);
    varid12 = netcdf.defVar(ncid,'Rad_770','double',[dimidx]);
    varid13 = netcdf.defVar(ncid,'Rad_800','double',[dimidx]);
    varid14 = netcdf.defVar(ncid,'Refl_680','double',[dimidx]);
    varid15 = netcdf.defVar(ncid,'Refl_705','double',[dimidx]);
    varid16 = netcdf.defVar(ncid,'Refl_740','double',[dimidx]);
    varid17 = netcdf.defVar(ncid,'Refl_770','double',[dimidx]);
    varid18 = netcdf.defVar(ncid,'Refl_800','double',[dimidx]);
    % Ha大气吸收窗口
    varid19 = netcdf.defVar(ncid,'SIF_Ha_sFLD','double',[dimidx]);
    varid20 = netcdf.defVar(ncid,'SIF_Ha_3FLD','double',[dimidx]);
    varid21 = netcdf.defVar(ncid,'SIF_Ha_iFLD','double',[dimidx]);
    varid22 = netcdf.defVar(ncid,'SIF_Ha_pFLD','double',[dimidx]);
    varid23 = netcdf.defVar(ncid,'SIF_Ha_SFM','double',[dimidx]);
    varid24 = netcdf.defVar(ncid,'SIF_Ha_SVD','double',[dimidx]);
    varid25 = netcdf.defVar(ncid,'SIF_Ha_DOAS','double',[dimidx]);
    % O2B大气吸收窗口
    varid26 = netcdf.defVar(ncid,'SIF_O2B_sFLD','double',[dimidx]);
    varid27 = netcdf.defVar(ncid,'SIF_O2B_3FLD','double',[dimidx]);
    varid28 = netcdf.defVar(ncid,'SIF_O2B_iFLD','double',[dimidx]);
    varid29 = netcdf.defVar(ncid,'SIF_O2B_pFLD','double',[dimidx]);
    varid30 = netcdf.defVar(ncid,'SIF_O2B_SFM','double',[dimidx]);
    varid31 = netcdf.defVar(ncid,'SIF_O2B_SVD','double',[dimidx]);
    varid32 = netcdf.defVar(ncid,'SIF_O2B_DOAS','double',[dimidx]);
    % H2O大气吸收窗口
    varid33 = netcdf.defVar(ncid,'SIF_H2O_sFLD','double',[dimidx]);
    varid34 = netcdf.defVar(ncid,'SIF_H2O_3FLD','double',[dimidx]);
    varid35 = netcdf.defVar(ncid,'SIF_H2O_iFLD','double',[dimidx]);
    varid36 = netcdf.defVar(ncid,'SIF_H2O_pFLD','double',[dimidx]);
    varid37 = netcdf.defVar(ncid,'SIF_H2O_SFM','double',[dimidx]);
    varid38 = netcdf.defVar(ncid,'SIF_H2O_SVD','double',[dimidx]);
    varid39 = netcdf.defVar(ncid,'SIF_H2O_DOAS','double',[dimidx]);
    % O2A大气吸收窗口
    varid40 = netcdf.defVar(ncid,'SIF_O2A_sFLD','double',[dimidx]);
    varid41 = netcdf.defVar(ncid,'SIF_O2A_3FLD','double',[dimidx]);
    varid42 = netcdf.defVar(ncid,'SIF_O2A_iFLD','double',[dimidx]);
    varid43 = netcdf.defVar(ncid,'SIF_O2A_pFLD','double',[dimidx]);
    varid44 = netcdf.defVar(ncid,'SIF_O2A_SFM','double',[dimidx]);
    varid45 = netcdf.defVar(ncid,'SIF_O2A_SVD','double',[dimidx]);
    varid46 = netcdf.defVar(ncid,'SIF_O2A_DOAS','double',[dimidx]);
    % 全波段反演结果
    varid47 = netcdf.defVar(ncid,'SIF_Fullband','double',[dimidx dimidy]);
    
    %---------------------------define attributes of the new variables--------------%  
    netcdf.putAtt(ncid,varid1,'units','/');                                                     %单位信息和long_name，其它的信息可依此定义
    netcdf.putAtt(ncid,varid2,'units','degree');      
    netcdf.putAtt(ncid,varid3,'units','nm');  
    netcdf.putAtt(ncid,varid4,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid5,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid6,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid7,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid8,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid9,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid10,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid11,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid12,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid13,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid14,'units','/'); 
    netcdf.putAtt(ncid,varid15,'units','/');
    netcdf.putAtt(ncid,varid16,'units','/');
    netcdf.putAtt(ncid,varid17,'units','/');
    netcdf.putAtt(ncid,varid18,'units','/');
    % Ha大气吸收窗口
    netcdf.putAtt(ncid,varid19,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid20,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid21,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid22,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid23,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid24,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid25,'units','mW/m2/nm/sr');
    % O2B大气吸收窗口
    netcdf.putAtt(ncid,varid26,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid27,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid28,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid29,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid30,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid31,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid32,'units','mW/m2/nm/sr');
    % H2O大气吸收窗口
    netcdf.putAtt(ncid,varid33,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid34,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid35,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid36,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid37,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid38,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid39,'units','mW/m2/nm/sr');
    % O2A大气吸收窗口
    netcdf.putAtt(ncid,varid40,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid41,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid42,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid43,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid44,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid45,'units','mW/m2/nm/sr');
    netcdf.putAtt(ncid,varid46,'units','mW/m2/nm/sr');
    % 全波段反演结果
    netcdf.putAtt(ncid,varid47,'units','mW/m2/nm/sr');
    
    netcdf.putAtt(ncid,varid1,'long_name','decimal Julian day');                                                     %单位信息和long_name，其它的信息可依此定义
    netcdf.putAtt(ncid,varid2,'long_name','Solar zenith angle');   
    netcdf.putAtt(ncid,varid3,'long_name','Wavelength');  
    netcdf.putAtt(ncid,varid4,'long_name','Solar irradiance / pi at 680 nm [4nm bin]');
    netcdf.putAtt(ncid,varid5,'long_name','Solar irradiance / pi at 705 nm [4nm bin]');
    netcdf.putAtt(ncid,varid6,'long_name','Solar irradiance / pi at 740 nm [4nm bin]');
    netcdf.putAtt(ncid,varid7,'long_name','Solar irradiance / pi at 770 nm [4nm bin]');
    netcdf.putAtt(ncid,varid8,'long_name','Solar irradiance / pi at 800 nm [4nm bin]');
    netcdf.putAtt(ncid,varid9,'long_name','Canopy radiance at 680 nm [4nm bin]');
    netcdf.putAtt(ncid,varid10,'long_name','Canopy radiance at 705 nm [4nm bin]');
    netcdf.putAtt(ncid,varid11,'long_name','Canopy radiance at 740 nm [4nm bin]');
    netcdf.putAtt(ncid,varid12,'long_name','Canopy radiance at 770 nm [4nm bin]');
    netcdf.putAtt(ncid,varid13,'long_name','Canopy radiance at 800 nm [4nm bin]');
    netcdf.putAtt(ncid,varid14,'long_name','Canopy reflectance at 680 nm [4nm bin]');
    netcdf.putAtt(ncid,varid15,'long_name','Canopy reflectance at 705 nm [4nm bin]');
    netcdf.putAtt(ncid,varid16,'long_name','Canopy reflectance at 740 nm [4nm bin]');
    netcdf.putAtt(ncid,varid17,'long_name','Canopy reflectance at 770 nm [4nm bin]');
    netcdf.putAtt(ncid,varid18,'long_name','Canopy reflectance at 800 nm [4nm bin]');
    % Ha大气吸收窗口
    netcdf.putAtt(ncid,varid19,'long_name','Canopy SIF at 657 nm by sFLD');
    netcdf.putAtt(ncid,varid20,'long_name','Canopy SIF at 657 nm by 3FLD');
    netcdf.putAtt(ncid,varid21,'long_name','Canopy SIF at 657 nm by iFLD');
    netcdf.putAtt(ncid,varid22,'long_name','Canopy SIF at 657 nm by pFLD');
    netcdf.putAtt(ncid,varid23,'long_name','Canopy SIF at 657 nm by SFM');
    netcdf.putAtt(ncid,varid24,'long_name','Canopy SIF at 657 nm by SVD');
    netcdf.putAtt(ncid,varid25,'long_name','Canopy SIF at 657 nm by DOAS');
    % O2B大气吸收窗口
    netcdf.putAtt(ncid,varid26,'long_name','Canopy SIF at 687 nm by sFLD');
    netcdf.putAtt(ncid,varid27,'long_name','Canopy SIF at 687 nm by 3FLD');
    netcdf.putAtt(ncid,varid28,'long_name','Canopy SIF at 687 nm by iFLD');
    netcdf.putAtt(ncid,varid29,'long_name','Canopy SIF at 687 nm by pFLD');
    netcdf.putAtt(ncid,varid30,'long_name','Canopy SIF at 687 nm by SFM');
    netcdf.putAtt(ncid,varid31,'long_name','Canopy SIF at 687 nm by SVD');
    netcdf.putAtt(ncid,varid32,'long_name','Canopy SIF at 687 nm by DOAS');
    % H2O大气吸收窗口
    netcdf.putAtt(ncid,varid33,'long_name','Canopy SIF at 719 nm by sFLD');
    netcdf.putAtt(ncid,varid34,'long_name','Canopy SIF at 719 nm by 3FLD');
    netcdf.putAtt(ncid,varid35,'long_name','Canopy SIF at 719 nm by iFLD');
    netcdf.putAtt(ncid,varid36,'long_name','Canopy SIF at 719 nm by pFLD');
    netcdf.putAtt(ncid,varid37,'long_name','Canopy SIF at 719 nm by SFM');
    netcdf.putAtt(ncid,varid38,'long_name','Canopy SIF at 719 nm by SVD');
    netcdf.putAtt(ncid,varid39,'long_name','Canopy SIF at 719 nm by DOAS');
    % O2A大气吸收窗口
    netcdf.putAtt(ncid,varid40,'long_name','Canopy SIF at 761 nm by sFLD');
    netcdf.putAtt(ncid,varid41,'long_name','Canopy SIF at 761 nm by 3FLD');
    netcdf.putAtt(ncid,varid42,'long_name','Canopy SIF at 761 nm by iFLD');
    netcdf.putAtt(ncid,varid43,'long_name','Canopy SIF at 761 nm by pFLD');
    netcdf.putAtt(ncid,varid44,'long_name','Canopy SIF at 761 nm by SFM');
    netcdf.putAtt(ncid,varid45,'long_name','Canopy SIF at 761 nm by SVD');
    netcdf.putAtt(ncid,varid46,'long_name','Canopy SIF at 761 nm by DOAS');
    % 全波段反演结果
    netcdf.putAtt(ncid,varid47,'long_name','Canopy SIF at full-wave band by F_SFM');
    
    netcdf.endDef(ncid);
    %
    netcdf.putVar(ncid,varid1,data.fracDOY);
    netcdf.putVar(ncid,varid2,data.SZA);
    netcdf.putVar(ncid,varid3,data.wl);
    netcdf.putVar(ncid,varid4,data.rad_sky_680);
    netcdf.putVar(ncid,varid5,data.rad_sky_705);
    netcdf.putVar(ncid,varid6,data.rad_sky_740);
    netcdf.putVar(ncid,varid7,data.rad_sky_770);
    netcdf.putVar(ncid,varid8,data.rad_sky_800);
    netcdf.putVar(ncid,varid9,data.rad_veg_680);
    netcdf.putVar(ncid,varid10,data.rad_veg_705);
    netcdf.putVar(ncid,varid11,data.rad_veg_740);
    netcdf.putVar(ncid,varid12,data.rad_veg_770);
    netcdf.putVar(ncid,varid13,data.rad_veg_800);
    netcdf.putVar(ncid,varid14,data.refl_680);
    netcdf.putVar(ncid,varid15,data.refl_705);
    netcdf.putVar(ncid,varid16,data.refl_740);
    netcdf.putVar(ncid,varid17,data.refl_770);
    netcdf.putVar(ncid,varid18,data.refl_800);
    % Ha大气吸收窗口
    netcdf.putVar(ncid,varid19,result.Ha.sFLD_SIF);
    netcdf.putVar(ncid,varid20,result.Ha.FLD3_SIF);
    netcdf.putVar(ncid,varid21,result.Ha.iFLD_SIF);
    netcdf.putVar(ncid,varid22,result.Ha.pFLD_SIF);
    netcdf.putVar(ncid,varid23,result.Ha.SFM_SIF);
    netcdf.putVar(ncid,varid24,result.Ha.SVD_SIF);
    netcdf.putVar(ncid,varid25,result.Ha.DOAS_SIF);
    % O2B大气吸收窗口
    netcdf.putVar(ncid,varid26,result.O2B.sFLD_SIF);
    netcdf.putVar(ncid,varid27,result.O2B.FLD3_SIF);
    netcdf.putVar(ncid,varid28,result.O2B.iFLD_SIF);
    netcdf.putVar(ncid,varid29,result.O2B.pFLD_SIF);
    netcdf.putVar(ncid,varid30,result.O2B.SFM_SIF);
    netcdf.putVar(ncid,varid31,result.O2B.SVD_SIF);
    netcdf.putVar(ncid,varid32,result.O2B.DOAS_SIF);
    % H2O大气吸收窗口
    netcdf.putVar(ncid,varid33,result.H2O.sFLD_SIF);
    netcdf.putVar(ncid,varid34,result.H2O.FLD3_SIF);
    netcdf.putVar(ncid,varid35,result.H2O.iFLD_SIF);
    netcdf.putVar(ncid,varid36,result.H2O.pFLD_SIF);
    netcdf.putVar(ncid,varid37,result.H2O.SFM_SIF);
    netcdf.putVar(ncid,varid38,result.H2O.SVD_SIF);
    netcdf.putVar(ncid,varid39,result.H2O.DOAS_SIF);
    % O2A大气吸收窗口
    netcdf.putVar(ncid,varid40,result.O2A.sFLD_SIF);
    netcdf.putVar(ncid,varid41,result.O2A.FLD3_SIF);
    netcdf.putVar(ncid,varid42,result.O2A.iFLD_SIF);
    netcdf.putVar(ncid,varid43,result.O2A.pFLD_SIF);
    netcdf.putVar(ncid,varid44,result.O2A.SFM_SIF);
    netcdf.putVar(ncid,varid45,result.O2A.SVD_SIF);
    netcdf.putVar(ncid,varid46,result.O2A.DOAS_SIF);
    % 全波段反演结果
    netcdf.putVar(ncid,varid47,result.FSFM.FSFM_SIF);
    
    netcdf.close(ncid);
end



