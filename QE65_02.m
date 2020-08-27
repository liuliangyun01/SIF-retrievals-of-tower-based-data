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

% ----------��Ҫ����: ����У��---------------------------------------------------------- %%
% --�˲��֣�ע���޸�
    Name_Tower = 'DM';                                    % ���޸�: station
    Year = '2018';                                        % ���޸�: year
% --·������
    Path_raw = strcat('.\Datasets\', Name_Tower, '\������\', Year, '\', Name_Tower, '_Uncor_Cor_Rad_Ref\');       % ���޸ģ�ԭʼ����·��
    Path_pro = strcat('.\Datasets\', Name_Tower, '\������\');                                                 % ���޸ģ�����������·��
    Path_Meteo =  strcat('.\Settings\', Name_Tower, '_', Year, '_��ѹ_�¶�_����.txt');                      % ���޸ģ���ѹ_�¶�����·��
    
switch Name_Tower
    case 'XTS'
        % --����������͸����·��
        Path_LUT_Dir = '.\Settings\Atmospheric_Transmittance\LLY_Ground_0050_SR_030_LUT_Dir_Tower_2.mat';      
        Path_LUT_Tot = '.\Settings\Atmospheric_Transmittance\LLY_Ground_0050_SR_030_LUT_Tot_Tower_2.mat';     
        % --��Ӧվ��ľ�γ��
        lon = 116.443;     lat = 40.179;
        TimeZone = 8;                     % ʱ����Ϊ8,��¼�ı���ʱ��
        % ����
        altitude = 50;
        % ����
        H_tower = 3.20;
        % ���˹۲����Ƕ�
        Angle = 25;
        % ��׼��ѹ���¶�
        P_zero = 1007.137;
        T_zero = 293.980;
    case 'HL'
        % --����������͸����·��
        Path_LUT_Dir = '.\Settings\Atmospheric_Transmittance\LLY_Ground_0500_SR_030_LUT_Dir_Tower_2.mat';      
        Path_LUT_Tot = '.\Settings\Atmospheric_Transmittance\LLY_Ground_0500_SR_030_LUT_Tot_Tower_2.mat';    
        % --��Ӧվ��ľ�γ��
        lon = 115.788;     lat = 40.349;
        TimeZone = 8;
        % ����
        altitude = 500;
        % ����
        H_tower = 4.08;
        % ���˹۲����Ƕ�
        Angle = 25;
        % ��׼��ѹ���¶�
        P_zero = 955.887;
        T_zero = 293.980;
    case 'DM'   
        % --����������͸����·��
        Path_LUT_Dir = '.\Settings\Atmospheric_Transmittance\LLY_Ground_1500_SR_030_LUT_Dir_Tower_20.mat';      
        Path_LUT_Tot = '.\Settings\Atmospheric_Transmittance\LLY_Ground_1500_SR_030_LUT_Tot_Tower_20.mat';    
        % --��Ӧվ��ľ�γ��
        lon = 100.372;     lat = 38.856;
        TimeZone = 8;
        % ����
        altitude = 1500;
        % ����
        H_tower = 24.5;
        % ���˹۲����Ƕ�
        Angle = 25;
        % ��׼��ѹ���¶�
        P_zero = 850.338;
        T_zero = 287.540;
        
    otherwise
        warning('Unexpected folder!')                                        % ���������ڴ�վ������
end

%% ------��ȡδУ�����������ݡ��Լ���ѹ�¶Ⱥ�͸��������------ %%
% --���ش���У��������ݣ���ȡ��ѹ���¶������Լ�������͸����
[Meteo] = Load_Meteo(Path_Meteo);                              % Load_Meteo�Ǻ�����������ѹ���¶����ݣ�Ϊ�˼����Ч����
[Table_up, Table_down] = Load_LUT(Path_LUT_Dir, Path_LUT_Tot); % Load_LUT�Ǻ��������ز��ұ�
% --����δУ���������ݣ������д���У��
files = dir(strcat(Path_raw, Name_Tower,'_Uncor_', '*.xlsx'));  % ����ֻ��ȡ��xlsx�ļ�
nfiles = size(files, 1);                              % xlsx�ļ�����


for i = 1:1:nfiles
    Path_i = strcat(Path_raw, files(i).name);       % �ļ�·��
    [~, name, ext] = fileparts(files(i).name);        % ��ȡָ���ļ���·�����ļ�������չ�����ļ���ʽ����
    Date = name(isstrprop(name, 'digit'));            % ��ȡ�ַ����е����֣����۲����ڣ�strrepҲ����ʵ��  
    
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
    
    % --����У������
    disp('Starting�� ����У��' + string(Date) + ', ֮�󱣴�');   tic;
    % --����ֲ���������߶ȣ������ڣ�
    [Time_init, Total_days, Veg_0, Veg_max] = Setup_Veg(Name_Tower, Year, Date);      % ע�⣺��Ӧվ���ֲ�����۲������������߶ȵ�
    veg_cor = [];
    sky_cor = [];
    for n = 1:1:data.nmeas
        % --ֲ���ĸ߶�____ע�⣺������ֲʱ�����������ո߶� 
        Year_i = str2num(strcat(Date(1), Date(2), Date(3), Date(4)));
        Month_i = str2num(strcat(Date(5), Date(6)));
        Day_i = str2num(strcat(Date(7), Date(8)));
        if Total_days == 0                              % ������жϣ��ܹ����ֿ�
            H_veg = 0;
        else
            % ����ֲ��ĳ��ʱ���ĸ߶�
            H_veg = ((datenum(Year_i, Month_i, Day_i) - Time_init) * (Veg_max-Veg_0) + Veg_0) ./ Total_days; 
        end
        
        % --������·�����ȣ� �۲�Ƕ�Angle��̫���춥��SZA
        L_up = ((H_tower - H_veg) ./ 1000) ./ cosd(Angle);                 
        L_down = ((H_tower - H_veg) ./ 1000) ./ cosd(data.sza(1, n));
        % --��Ч·��
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
        % --E790/E660, 660��790�ֱ��ʾ��660.044nm��790.113nm
        Index790 = data.wl == 790.113; 
        Index660 = data.wl == 660.044; 
        E_Ratio = data.sky(Index790, n) ./ data.sky(Index660, n);
        % --���������е�͸����, LUT_Tra����
        [Up_tra] = LUT_Tra(Table_up, E_Ratio, RTPL_up);               % LUT_Tra�Ǻ���������͸����
        [Down_tra] = LUT_Tra(Table_down, E_Ratio, RTPL_down);
        New_Up_tra = [Up_tra(1,1:20), Up_tra,Up_tra(1,(end-20+1):end)];
        New_Down_tra = [Down_tra(1,1:20), Down_tra, Down_tra(1,(end-20+1):end)];
        [~, Up_tra_min_num] = min(New_Up_tra(1,:), [], 2);
        [~, ~] = min(New_Down_tra(1,:), [], 2);
        Num_tra = Up_tra_min_num;                                          % ����͸���ʵ���Сֵ���ȶ�
        Window_O2A = [758, 770]; 
        [Range, Index_veg_min, Index_sky_min] = Cal_Index(Window_O2A, data); % Cal_IndexΪ�������ҵ�������������ߣ�����A������Сֵ�������
        Num_qe = (Index_veg_min + Index_sky_min) ./ 2;
        Inter_bef = Num_tra - (Num_qe - 1);
        Inter_aft = Num_tra + (1021 - Num_qe);
        
        % --��TOA(��������)��ֵת��ΪTOC(�ڲ㴦)��ֵ-----
        veg_cor(:, n) = data.veg(:,n) ./ New_Up_tra(1,Inter_bef:Inter_aft)';
        sky_cor(:, n) = data.sky(:,n) .* New_Down_tra(1,Inter_bef:Inter_aft)';
    end
    
    data.veg_cor = veg_cor;
    data.sky_cor = sky_cor;
    data.ref_cor = data.veg_cor ./ data.sky_cor;

    % --��������Ϊnc, xlsx��mat��ʽ
    save_xlsx(Path_pro, Name_Tower, Date, data);
    save_mat(Path_pro, Name_Tower, Date, data);
    save_nc(Path_pro, Name_Tower, Date, data);
    
end
disp('Time delays: ' + string(toc) + 's');


%% ------------------------��������------------------------ %%
% --������ע���жϴ��������պ͵��������Сֵ�����к�
function [Range, Index_veg_min_ave, Index_sky_min_ave] = Cal_Index(wl_range, data)
    Range = find(data.wl > wl_range(1) & data.wl < wl_range(2));    % ���ղ��η�Χ������
    Index_left = Range(1);
    Index_right = Range(end);
    for i = 1:1:data.nmeas
        [veg_min(1,i), Index_veg_min(1,i)] = min(data.veg(Index_left:Index_right, i),[],1);
        [sky_min(1,i), Index_sky_min(1,i)] = min(data.sky(Index_left:Index_right, i),[],1);
    end
    % ���ա��������Сֵ�ڳ�������������ݴ����쳣����ȥ����������10%����������� %
    Bound_min = round(data.nmeas*0.1);
    Bound_max = round(data.nmeas)-round(data.nmeas*0.1);
    Index_veg_min_ave = round(mean(Index_veg_min(1,Bound_min:Bound_max),2))+ Index_left-1;
    Index_sky_min_ave = round(mean(Index_sky_min(1,Bound_min:Bound_max),2))+ Index_left-1;
end

% --������ע��������ѹ���¶����ݣ�����ʽת�����磺2017/5/1 0:20ת��Ϊ��20170501 0.01389��
function [Meteo] = Load_Meteo(Path_Meteo)

    [Date, Time, Tem, Press] = textread(Path_Meteo, '%s%s%s%s', 'headerlines', 1);
    % ----------�����+���ڡ���ʽת��---------- %
    nyear = size(Date, 1);
    for i = 1:1:nyear
        Date_i = Date{i,1};
        Length_Date_i = length(Date_i);
        Num_dia = strfind(Date_i,'/');
        % ��ȡ������
        Year = strcat(Date_i(1), Date_i(2), Date_i(3), Date_i(4));
        % ��ȡ������
        if Num_dia(1,2) - Num_dia(1,1) == 2
            Month = strcat('0', Date_i(Num_dia(1,1)+1));
        elseif Num_dia(1,2) - Num_dia(1,1) == 3
            Month = strcat(Date_i(Num_dia(1,1)+1), Date_i(Num_dia(1,1)+2));
        end
        % ��ȡ������
        if Length_Date_i - Num_dia(1,2) == 1
            Day = strcat('0', Date_i(Num_dia(1,2)+1));
        elseif Length_Date_i - Num_dia(1,2) == 2
            Day = strcat(Date_i(Num_dia(1,2)+1), Date_i(Num_dia(1,2)+2));
        end
        Date_new(i,1) = str2double(strcat(Year, Month, Day));
    end
    
    % ----------��ʱ+�֡���ʽת��---------- %
    ntime = size(Time, 1);
    for j = 1:1:ntime
        Time_j = Time{j, 1};
        Num_col = strfind(Time_j, ':');
        %��ȡʱ��������
        if Num_col(1,1) == 2
            Hour = str2double(Time_j(Num_col(1,1)-1));
            Minute = str2double(strcat(Time_j(Num_col(1,1)+1), Time_j(Num_col(1,1)+2)));
        elseif Num_col(1,1) == 3
            Hour = str2double(strcat(Time_j(Num_col(1,1)-2), Time_j(Num_col(1,1)-1)));
            Minute = str2double(strcat(Time_j(Num_col(1,1)+1), Time_j(Num_col(1,1)+2)));
        end
        Time_new(j,1) = Hour./24 + Minute./60./24;
    end
    
    % ----------�¶ȡ���ѹ���ַ���ʽת��Ϊ���ָ�ʽ---------- %
    ntem = size(Tem, 1);
    for k = 1:1:ntem
        Tem_new(k, 1) = str2double(Tem{k,1});
        Press_new(k, 1) = str2double(Press{k,1});
    end

    Meteo = [Date_new, Time_new, Tem_new, Press_new];
end

% --������ע������������͸���ʣ�MORTRANģ������
function [Table_up_new, Table_down_new] = Load_LUT(Path_LUT_Dir, Path_LUT_Tot)
    disp('Starting�� ����������͸���� ');
    % ----��ȡE��ֵ _Tdir��ģ������ && _Ttot��ģ������---- %
    LUT_Dir = load(Path_LUT_Dir);      Name_Dir = whos('-file', Path_LUT_Dir);  % ����֪���ṹ������������
    LUT_Tot = load(Path_LUT_Tot);      Name_Tot = whos('-file', Path_LUT_Tot); 
    Table_up_raw = LUT_Dir.(Name_Dir.name);
    Table_down_raw = LUT_Tot.(Name_Tot.name); 
    AOD_00_90_up = Table_up_raw(:,2);  
    AOD_00_90_down = Table_down_raw(:,2);
    for i = 1:1:10             % ������\���� ģ�����ݰ���AOD�̶���VZA��С��������
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

% ---������ע��ͨ�����ұ�ķ�ʽ������õ�����͸����
function [Transmittance] = LUT_Tra(Table_up_new, E_Ratio, RTPL_up)

    % ---��H_tower/cos(VZA)ȷ������E��ֵ��Tdir͸���ʵ���ֵ�������Բ�ֵ
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

    % ---��E��ֵȷ������Tdir͸���ʵ���ֵ�������Բ�ֵ
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

% ---������ע����������
function save_xlsx(Path_pro, Name_Tower, Date, data)
    disp('Saving�� xlsx');  
    
    % --------------�����ļ���-------------- %
    folder_date = strcat(Path_pro, Date(1:4), '\');    
    if exist(folder_date, 'dir')==0     % �жϡ���ݡ��ļ����Ƿ����
        mkdir(folder_date);             % ������ʱ�򣬴����ļ���
    end
    folder = strcat(folder_date, Name_Tower, '_Uncor_Cor_Rad_Ref\');     
    if exist(folder, 'dir')==0          % �жϡ�_Uncor_Rad_Ref_���ļ����Ƿ����
        mkdir(folder);                  % ������ʱ�򣬴����ļ���
    end
    % --------------���Excel���-------------- %
    Output_path = strcat(folder, Name_Tower, '_Cor_Rad_Ref_', Date, '.xlsx');       % �޸�:���.xlsx�ļ��еĹ۲�վ���ƣ��磺DM
    Out_Ref_time_sza = [data.time; data.sza];
    Out_Rad_time_sza = [data.time, data.time; data.sza, data.sza];
    Out_rad_veg_sky = [data.veg_cor, data.sky_cor];
    Names = {'ʱ��'; 'SZA'};
    % --���������յķ����ȣ�δУ����-- %
    xlswrite(Output_path, Names, '���������յķ����ȣ���У����', 'A1:A2');
    xlswrite(Output_path, data.wl, '���������յķ����ȣ���У����', 'A3');
    xlswrite(Output_path, Out_Rad_time_sza, '���������յķ����ȣ���У����', 'B1');
    xlswrite(Output_path, Out_rad_veg_sky, '���������յķ����ȣ���У����', 'B3');
    % --�����ʣ�δУ����-- %
    xlswrite(Output_path, Names, '�����ʣ���У����', 'A1:A2');
    xlswrite(Output_path, data.wl, '�����ʣ���У����', 'A3');
    xlswrite(Output_path, Out_Ref_time_sza, '�����ʣ���У����', 'B1');
    xlswrite(Output_path, data.ref_cor, '�����ʣ���У����', 'B3');
    
end
function save_mat(Path_pro, Name_Tower, Date, data)
    disp('Saving�� mat'); 

    % --------------�����ļ���-------------- %
    folder_date = strcat(Path_pro, Date(1:4), '\');    
    if exist(folder_date, 'dir')==0     
        mkdir(folder_date);            
    end
    folder = strcat(folder_date, Name_Tower, '_Uncor_Cor_Rad_Ref\');  
    if exist(folder, 'dir')==0    
        mkdir(folder);             
    end
    % --------------���mat���-------------- %
    Output_path = strcat(folder, Name_Tower, '_Cor_Rad_Ref_', Date, '.mat');       % ���
    
    wl = data.wl;
    time = data.time;
    sza = data.sza;
    veg = data.veg_cor;
    sky = data.sky_cor;
    ref = data.ref_cor;
    save(Output_path, 'wl', 'time', 'sza', 'veg', 'sky', 'ref');
end
function save_nc(Path_pro, Name_Tower, Date, data)
    disp('Saving�� nc'); 

    % --------------�����ļ���-------------- %
    folder_date = strcat(Path_pro, Date(1:4), '\');    
    if exist(folder_date, 'dir')==0     % �жϡ���ݡ��ļ����Ƿ����
        mkdir(folder_date);             % ������ʱ�򣬴����ļ���
    end
    folder = strcat(folder_date, Name_Tower, '_Uncor_Cor_Rad_Ref\');     
    if exist(folder, 'dir')==0          % �жϡ�_Uncor_Rad_Ref_���ļ����Ƿ����
        mkdir(folder);                  % ������ʱ�򣬴����ļ���
    end
    % --------------���nc���-------------- %
    Output_path = strcat(folder, Name_Tower, '_Cor_Rad_Ref_', Date, '.nc');       % �޸�:���.xlsx�ļ��еĹ۲�վ���ƣ��磺DM
    
    % ----------DEFINE THE FILE---------- %                                        
    ncid = netcdf.create(Output_path,'CLOBBER');   % ����һ��������ݵ�nc�ļ�
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
    netcdf.putAtt(ncid,varid1,'units','/');       %��λ��Ϣ��long_name����������Ϣ�����˶���
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
    % �������
    netcdf.putVar(ncid,varid1,data.time);
    netcdf.putVar(ncid,varid2,data.sza);
    netcdf.putVar(ncid,varid3,data.wl);
    netcdf.putVar(ncid,varid4,data.veg_cor);
    netcdf.putVar(ncid,varid5,data.sky_cor);
    netcdf.putVar(ncid,varid6,data.ref_cor);
    
    netcdf.close(ncid);
end
