function [windows] = Setup_Window()
    
    windows = [];
    % -----Attention: fitting window range for sFLD、3FLD、SFM、SVD、DOAS;----- %
    % 用于sFLD, iFLD, pFLD
    FLD = [];
    FLD.outerH = [656.0, 656.2];     % 端点波段，取均值
    FLD.outerB = [687.0, 687.2];   
    FLD.outerW = [717.0, 717.2];  
    FLD.outerA = [758.6, 758.8]; 
    FLD.windowH = [656, 659];         % 范围内求最小值
    FLD.windowB = [686, 689];    
    FLD.windowW = [717, 722]; 
    FLD.windowA = [758, 770]; 
    windows.FLD = FLD;

    % 3FLD
    FLD3 = [];
    FLD3.FLD3_leftH  = [656.0, 656.2];    FLD3.FLD3_rightH = [659.0, 659.2];  
    FLD3.FLD3_leftB  = [687.0, 687.2];    FLD3.FLD3_rightB = [689.0, 689.2];    
    FLD3.FLD3_leftW  = [717.0, 717.2];    FLD3.FLD3_rightW = [722.0, 722.2]; 
    FLD3.FLD3_leftA  = [758.6, 758.8];    FLD3.FLD3_rightA = [770.0, 770.2];
    FLD3.windowH = [656, 659];           % 范围内求最小值
    FLD3.windowB = [686, 689];    
    FLD3.windowW = [717, 722]; 
    FLD3.windowA = [758, 770]; 
    windows.FLD3 = FLD3;
    
    % SFM
    SFM = [];
    SFM.SFM_windowH = [656, 667];       % 波段范围
    SFM.SFM_windowB = [687, 698];
    SFM.SFM_windowW = [717, 728];
    SFM.SFM_windowA = [759, 770];
    SFM.SFM_Ha =  [656.0 656.2];       % 某处的荧光值,SFM和SVD算法直接定死了
    SFM.SFM_O2B = [687.0 687.2];
    SFM.SFM_H2O = [720.0 720.2];
    SFM.SFM_O2A = [760.0 760.2];
    windows.SFM = SFM;
    % SVD
    SVD = [];
    SVD.SVD_windowH  = [655, 670];     
    SVD.SVD_windowB  = [675, 710];      % wide O2B and O2A fitting window; % e.g. narrow O2A fitting window:[759.86 762.79];
    SVD.SVD_windowW  = [700, 740];
    SVD.SVD_windowA  = [740, 780];  
    SVD.SVD_Ha =  [656.0 656.2];       % 某处的荧光值,SFM和SVD算法直接定死了
    SVD.SVD_O2B = [687.0 687.2];
    SVD.SVD_H2O = [720.0 720.2];
    SVD.SVD_O2A = [760.0 760.2];
    % Example SVD settings
    % if rolling > 0, abs(rolling) = the number of spectra in a centered moving window used for training the SVD
    % if rolling = 0, it will use all spectra from the day for training the SVD
    % if rolling < 0, it will use that number of spectra but evenly dispersed throughout the day for training the SVD
    SVD.rolling = 0;                    % moving window of 5 spectra 
    windows.SVD = SVD;
    
    % DOAS, F_SFM, 主成分分析，也用于取波段范围内辐亮度最小值       
    Others = [];                           % 四种方法，暂用的波段范围，取最小值及序号 
    Others.WindowH = [656, 659];           % Attention: 注意修改
    Others.WindowB = [686, 689];    
    Others.WindowW = [717, 722]; 
    Others.WindowA = [758, 770]; 
    windows.Others  = Others;
    
    % -----Attention: Modify the range of spectrum;----- %
    % --得到下面五个，辐照度，辐亮度和反射率
    windows.wl_680 = [676, 680];    % 676-680nm, red
    windows.wl_705 = [703, 707];    % 703-707nm, red
    windows.wl_740 = [736, 740];    % 736-740nm, rededge
    windows.wl_770 = [770, 774];    % 770-774nm, near_red
    windows.wl_800 = [795, 799];    % 795-799nm, near_red
end