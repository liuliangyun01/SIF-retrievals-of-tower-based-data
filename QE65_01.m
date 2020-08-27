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

% ----------��Ҫ����: organize the observed radiance data in txt format into xlsx and mat format---------- %%
% --�˲��֣�ע���޸�
    Name_Tower = 'DM';                                    % ���޸�: station
    Year = '2018';                                        % ���޸�: year
% --·������
    Path_raw = strcat('.\Datasets\', Name_Tower, '\ԭʼ����\', Year, '\');    % ���޸ģ�ԭʼ����·�������Ƚ������ļ���, ��ԭʼ���ݷ������У�
    Path_pro = strcat('.\Datasets\', Name_Tower, '\������\');               % ���޸ģ�����������·�������Ƚ������ļ��У�
    Norm_wl(:, 1) = importdata('.\Settings\wl_pro.txt');                 % ���޸ģ���׼����·��
    
switch Name_Tower
    case 'XTS'
        % --��Ӧվ����������
        QEPB_serial = 'QEPB1220';
        if strcmp(Name_Tower, 'XTS') && strcmp(Year,'2017')                      %  --***ע�⣺2017��С��ɽ�ͻ�������,7��6��---- %                                             
            QEPB_serial = 'QEPB1218';
        end
        % --��Ӧվ��ľ�γ��
        lon = 116.443;     lat = 40.179;
        TimeZone = 8;                     % ʱ����Ϊ8,��¼�ı���ʱ��
    case 'HL'     
        % --��Ӧ���������
        QEPB_serial = 'QEPB1218';
        if strcmp(Name_Tower, 'HL') && strcmp(Year,'2017')                        %  --***ע�⣺2017��С��ɽ�ͻ�������---- %                                             
            QEPB_serial = 'QEPB1220';
        elseif strcmp(Name_Tower, 'HL') && strcmp(Year,'2020')                    %  --***ע�⣺2020�갢��ͻ�����������---- %                                             
            QEPB_serial = 'QEPB1564';
        end
        % --��Ӧվ��ľ�γ��
        lon = 115.788;     lat = 40.349;
        TimeZone = 8;
    case 'DM'            
        % --��Ӧ���������
        QEPB_serial = 'QEPB1219';
        % --��Ӧվ��ľ�γ��
        lon = 100.372;     lat = 38.856;
        TimeZone = 8;

    otherwise
        warning('Unexpected folder!')                                        % ���������ڴ�վ������
end

%% ------���ع۲����_txt��ʽ�ı��ļ����Ա�׼���β�ֵ------ %%
Files = dir(fullfile(Path_raw, strcat(Year, '*')));  % Struct of files
Length = size(Files, 1);                             % Number of daily observation data,days

for i = 1:1:Length
    
    date_i = Files(i,1).name;
    year_i = str2double(strcat(date_i(1), date_i(2), date_i(3), date_i(4)));
    month_i = str2double(strcat(date_i(5), date_i(6)));
    day_i = str2double(strcat(date_i(7), date_i(8)));
    if strcmp(Name_Tower, 'XTS') && (str2num(date_i) >= 20170706)  % ��������
        QEPB_serial = 'QEPB1220';
    elseif strcmp(Name_Tower, 'HL') && (str2num(date_i) >= 20170711)  % ��������
        QEPB_serial = 'QEPB1218';
    end
    path_i = strcat(Path_raw, Files(i,1).name, '\Auto\', QEPB_serial, '\'); % �����ļ���·��
    disp('Starting�� ��ȡ' + string(date_i) + '��������txt�ļ������Ա�׼���β�ֵ, ֮�󱣴�');   tic;
    %--��ȡ�����ļ��е�txt����
    txts = dir(fullfile(path_i, '*.txt'));     % �������е�txt�ļ����ṹ����ʽ
    ntxts = size(txts, 1);                     % ���죬txt������   
    wl = [];                                   % **----�������Ҫ���������ݣ���Ȼ�����ڴ�
    time = [];
    sza = [];
    spec = [];
    for j = 1:1:ntxts
        path_j = strcat(path_i, txts(j, 1).name);                           % ����txt��ʽ�ļ�·��
        [wl_j, spec_j, time_j] = get_spec(path_j);
        sza_j = Cal_SZA(time_j, year_i, month_i, day_i, lon, lat, TimeZone);
        wl(:,j) = wl_j(:,1);
        sza(:,j) = sza_j(:,1);
        time(:,j) = time_j(:,1);
        spec(:,j) = spec_j(:,1);
    end
    % --�����ݰ���ʱ����������-------- %
    [time, index] = sort(time);
    sza = sza(index);
    spec = spec(:, index);
    % --�Ե��������������-------- %
    time_mean = [];
    sza_mean = [];
    veg = [];
    sky = [];
    % --�ж�txt���ݵ������Ƿ�Ϊ3�ı���-------- %
    if mod(ntxts, 3) == 0                   
        for k = 1:1:(ntxts/3)
            time_mean(:, k) = mean(time(1, (3*k-2):3*k), 2);
            sza_mean(:, k) = mean(sza(1, (3*k-2):3*k), 2);
            veg(:, k) = spec(:, 3*k-1);
            % --�޳����������������-- %,    --{�����������ݣ����ڸ���}--
            % ***ע���޸�***���������Ч�����á�Irr_760/Irr_758��
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
    % ---�Ա�׼���β�ֵ---------------%
    Prim_wl(:,1) = wl(:,1);
    veg_interp = [];
    sky_interp = [];
    for m = 1:1:(ntxts/3)  
        veg_interp(:,m) = interp1(Prim_wl(:, 1), veg(:, m), Norm_wl(:,1), 'spline');
        sky_interp(:,m) = interp1(Prim_wl(:, 1), sky(:, m), Norm_wl(:,1), 'spline');
    end
    % **ע������___����Ҵ��20170511-20170711Ϊ�������ݣ���Ҫת����λ�������ա�������ֵ��λת����̫�������ȡ���������ȡ�----%
    if (str2double(date_i) >= 20170511) && (str2double(date_i) <= 20170711)
        veg_interp = veg_interp*10;
        sky_interp = sky_interp*10/pi;    % ע��***�����ﶨ��ϵ�������⣨by����濣�
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
    
    % --��������Ϊnc, xlsx��mat��ʽ
    save_xlsx(Path_pro, Name_Tower, date_i, data);
    save_mat(Path_pro, Name_Tower, date_i, data);
    save_nc(Path_pro, Name_Tower, date_i, data);
    
    disp('Time delays: ' + string(toc) + 's');
end

%% ------------------------��������------------------------ %%
% --������ע����ȡ��������
function [wl, spec, time] = get_spec(path)
    % --��ȡʱ��
    if nargin == 1                                         % �ж���������ĸ�����
        fid = fopen(path);                                 % fidΪfopen����ص��ļ���ʶ��
        a = textscan(fid, '%s%s%s%s', 'headerlines',10);   % �ҵ�starting
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
        time = (Hour/24.0 + Min/24.0/60.0 + Sec/24.0/60.0/60.0);   % ��ʱ�任���0-1��ʱ�䷶Χ
    end
    
    % --��ȡ����
    fid = fopen(path);
    b = textscan(fid, '%s%s%s%s', 'headerlines', 18);
    fclose(fid);
    b1 = b{1, 1};
    b4 = b{1, 4};
    N = size(b1,1);
    
    wl = str2num(char(b1));    % str2num,��str2doubleֻ�����ڱ�����ת��spec=str2double(b1);
    spec = str2num(char(b4));     
end

% --������ע����������Ϊxlsx��mat
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
    Output_path = strcat(folder, Name_Tower, '_Uncor_Rad_Ref_', Date, '.xlsx');       % �޸�:���.xlsx�ļ��еĹ۲�վ���ƣ��磺DM
    
    Out_Ref_time_sza = [data.time; data.sza];
    Out_Rad_time_sza = [data.time, data.time; data.sza, data.sza];
    Out_rad_veg_sky  = [data.veg, data.sky];
    Names = {'ʱ��'; 'SZA'};
    % --���������յķ����ȣ�δУ����-- %
    xlswrite(Output_path, Names, '���������յķ����ȣ�δУ����', 'A1:A2');
    xlswrite(Output_path, data.wl, '���������յķ����ȣ�δУ����', 'A3');
    xlswrite(Output_path, Out_Rad_time_sza, '���������յķ����ȣ�δУ����', 'B1');
    xlswrite(Output_path, Out_rad_veg_sky, '���������յķ����ȣ�δУ����', 'B3');

    % --�����ʣ�δУ����-- %
    xlswrite(Output_path, Names, '�����ʣ�δУ����', 'A1:A2');
    xlswrite(Output_path, data.wl, '�����ʣ�δУ����', 'A3');
    xlswrite(Output_path, Out_Ref_time_sza, '�����ʣ�δУ����', 'B1');
    xlswrite(Output_path, data.ref, '�����ʣ�δУ����', 'B3');

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
    Output_path = strcat(folder, Name_Tower, '_Uncor_Rad_Ref_', Date, '.mat');       % ���
    
    wl = data.wl;
    time = data.time;
    sza = data.sza;
    veg = data.veg;
    sky = data.sky;
    ref = data.ref;
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
    Output_path = strcat(folder, Name_Tower, '_Uncor_Rad_Ref_', Date, '.nc');       % �޸�:���.xlsx�ļ��еĹ۲�վ���ƣ��磺DM
    
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
    netcdf.putAtt(ncid,varid4,'long_name','Canopy radiance');
    netcdf.putAtt(ncid,varid5,'long_name','Solar irradiance / pi');
    netcdf.putAtt(ncid,varid6,'long_name','Canopy reflectance');

    netcdf.endDef(ncid);
    % �������
    netcdf.putVar(ncid,varid1,data.time);
    netcdf.putVar(ncid,varid2,data.sza);
    netcdf.putVar(ncid,varid3,data.wl);
    netcdf.putVar(ncid,varid4,data.veg);
    netcdf.putVar(ncid,varid5,data.sky);
    netcdf.putVar(ncid,varid6,data.ref);
    
    netcdf.close(ncid);
end

