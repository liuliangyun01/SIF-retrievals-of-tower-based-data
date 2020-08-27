% SIF retrieval wrapper
function [result] = SIFretrieval(data, models, Scope, windows)
    
    result = [];
    Ha = [];  O2B = [];   H2O = [];  O2A = [];  FSFM = [];
    % Flag，去除吸收波段的值，用于插值，或者主成分分析，或者全波段的反演
    Flags = [];
    Flags.Flag_H  = data.wl > windows.Others.WindowH(1) & data.wl < windows.Others.WindowH(2);
    Flags.Flag_B  = data.wl > windows.Others.WindowB(1) & data.wl < windows.Others.WindowB(2);
    Flags.Flag_W  = data.wl > windows.Others.WindowW(1) & data.wl < windows.Others.WindowW(2);
    Flags.Flag_A  = data.wl > windows.Others.WindowA(1) & data.wl < windows.Others.WindowA(2);
    Flags.Flag = Flags.Flag_H + Flags.Flag_B + Flags.Flag_W + Flags.Flag_A;
    
    % -----------------------各种算法------------------------ %
    % sFLD
    if models.sFLD == 1
        disp("Starting sFLD_"+ string(unique(floor(data.fracDOY))));                              % disp, output data
        tic;                                                % tic is used to save the current time, and then toc is used to record the program completion time.
        for i = 1:1:4
            Name_i = cell2mat(fieldnames(windows.Others));
            if strcmp(Name_i(i, end), 'H')   
                outer = windows.FLD.outerH;
                wl_range = windows.FLD.windowH;
                [sFLD_SIF, sFLD_ref] = sFLD(data, outer, wl_range);
                Ha.sFLD_SIF = sFLD_SIF;
                Ha.sFLD_ref = sFLD_ref;
            elseif strcmp(Name_i(i, end), 'B')   
                outer = windows.FLD.outerB;
                wl_range = windows.FLD.windowB;
                [sFLD_SIF, sFLD_ref] = sFLD(data, outer, wl_range);
                O2B.sFLD_SIF = sFLD_SIF;
                O2B.sFLD_ref = sFLD_ref;
            elseif strcmp(Name_i(i, end), 'W')  
                outer = windows.FLD.outerW;
                wl_range = windows.FLD.windowW;
                [sFLD_SIF, sFLD_ref] = sFLD(data, outer, wl_range);
                H2O.sFLD_SIF = sFLD_SIF;
                H2O.sFLD_ref = sFLD_ref;
            elseif strcmp(Name_i(i, end), 'A')   
                outer = windows.FLD.outerA;
                wl_range = windows.FLD.windowA;
                [sFLD_SIF, sFLD_ref] = sFLD(data, outer, wl_range);
                O2A.sFLD_SIF = sFLD_SIF;
                O2A.sFLD_ref = sFLD_ref;
            end
            
        end
        msg = 'Time delays for sFLD: ' + string(toc) + 's';   disp(msg);
    end
    % 3FLD
    if models.FLD3 == 1
        disp("Starting 3FLD_"+ string(unique(floor(data.fracDOY))));      tic;
        for i = 1:1:4
            Name_i = cell2mat(fieldnames(windows.Others));
            if strcmp(Name_i(i, end), 'H')   
                left = windows.FLD3.FLD3_leftH;     % 波谷左端
                right = windows.FLD3.FLD3_rightH;   % 波谷右端
                wl_range = windows.FLD3.windowH;
                [FLD3_SIF, iFLD_ref] = FLD3(data, left, right, wl_range);   
                Ha.FLD3_SIF = FLD3_SIF;
                Ha.FLD3_ref = iFLD_ref;
            elseif strcmp(Name_i(i, end), 'B')   
                left = windows.FLD3.FLD3_leftB;     % 波谷左端
                right = windows.FLD3.FLD3_rightB;   % 波谷右端
                wl_range = windows.FLD3.windowB;
                [FLD3_SIF, iFLD_ref] = FLD3(data, left, right, wl_range);   
                O2B.FLD3_SIF = FLD3_SIF;
                O2B.FLD3_ref = iFLD_ref;
            elseif strcmp(Name_i(i, end), 'W')  
                left = windows.FLD3.FLD3_leftW;     % 波谷左端
                right = windows.FLD3.FLD3_rightW;   % 波谷右端
                wl_range = windows.FLD3.windowW;
                [FLD3_SIF, iFLD_ref] = FLD3(data, left, right, wl_range);   
                H2O.FLD3_SIF = FLD3_SIF;
                H2O.FLD3_ref = iFLD_ref;
            elseif strcmp(Name_i(i, end), 'A')   
                left = windows.FLD3.FLD3_leftA;     % 波谷左端
                right = windows.FLD3.FLD3_rightA;   % 波谷右端
                wl_range = windows.FLD3.windowA;
                [FLD3_SIF, iFLD_ref] = FLD3(data, left, right, wl_range);   
                O2A.FLD3_SIF = FLD3_SIF;
                O2A.FLD3_ref = iFLD_ref;
            end
        end
        
        msg = 'Time delays for 3FLD: ' + string(toc) + 's';   disp(msg);
    end
    % iFLD & pFLD
    if models.iFLD == 1 && models.pFLD == 1
        disp("Starting iFLD & pFLD_"+ string(unique(floor(data.fracDOY))));      tic;
        for i = 1:1:4
            Name_i = cell2mat(fieldnames(windows.Others));
            if strcmp(Name_i(i, end), 'H')   
                outer = windows.FLD.outerH;
                wl_range = windows.FLD.windowH;
                [iFLD_SIF, iFLD_ref] = iFLD(data, outer, wl_range, Flags);   
                Ha.iFLD_SIF = iFLD_SIF;
                Ha.iFLD_ref = iFLD_ref;
                [pFLD_SIF, pFLD_ref] = pFLD(data, outer, wl_range, Flags, Scope); 
                Ha.pFLD_SIF = pFLD_SIF;
                Ha.pFLD_ref = pFLD_ref;
            elseif strcmp(Name_i(i, end), 'B')   
                outer = windows.FLD.outerB;
                wl_range = windows.FLD.windowB;
                [iFLD_SIF, iFLD_ref] = iFLD(data, outer, wl_range, Flags);   
                O2B.iFLD_SIF = iFLD_SIF;
                O2B.iFLD_ref = iFLD_ref;
                [pFLD_SIF, pFLD_ref] = pFLD(data, outer, wl_range, Flags, Scope); 
                O2B.pFLD_SIF = pFLD_SIF;
                O2B.pFLD_ref = pFLD_ref;
            elseif strcmp(Name_i(i, end), 'W')  
                outer = windows.FLD.outerW;
                wl_range = windows.FLD.windowW;
                [iFLD_SIF, iFLD_ref] = iFLD(data, outer, wl_range, Flags);   
                H2O.iFLD_SIF = iFLD_SIF;
                H2O.iFLD_ref = iFLD_ref;
                [pFLD_SIF, pFLD_ref] = pFLD(data, outer, wl_range, Flags, Scope); 
                H2O.pFLD_SIF = pFLD_SIF;
                H2O.pFLD_ref = pFLD_ref;
            elseif strcmp(Name_i(i, end), 'A')   
                outer = windows.FLD.outerA;
                wl_range = windows.FLD.windowA;
                [iFLD_SIF, iFLD_ref] = iFLD(data, outer, wl_range, Flags);   
                O2A.iFLD_SIF = iFLD_SIF;
                O2A.iFLD_ref = iFLD_ref;
                [pFLD_SIF, pFLD_ref] = pFLD(data, outer, wl_range, Flags, Scope); 
                O2A.pFLD_SIF = pFLD_SIF;
                O2A.pFLD_ref = pFLD_ref;
            end
        end
        msg = 'Time delays for iFLD & pFLD: ' + string(toc) + 's';   disp(msg);
    end
    % SFM ,assume polynomial refl and fluo, 2nd order poly  
    if models.SFM == 1 
        disp("Starting SFM_"+ string(unique(floor(data.fracDOY))));      tic;
        for i = 1:1:4
            Name_i = cell2mat(fieldnames(windows.Others));
            if strcmp(Name_i(i, end), 'H')   
                sfm_wlrange = windows.SFM.SFM_windowH;
                sfm_Ha = windows.SFM.SFM_Ha;
                [fluo, sfm] = SFM(data, sfm_wlrange, sfm_Ha);   
                Ha.SFM_SIF = fluo;
                Ha.SFM = sfm;
            elseif strcmp(Name_i(i, end), 'B')   
                sfm_wlrange = windows.SFM.SFM_windowB;
                sfm_O2B = windows.SFM.SFM_O2B;
                [fluo, sfm] = SFM(data, sfm_wlrange, sfm_O2B);   
                O2B.SFM_SIF = fluo;
                O2B.SFM = sfm;
            elseif strcmp(Name_i(i, end), 'W')   
                sfm_wlrange = windows.SFM.SFM_windowW;
                sfm_H2O = windows.SFM.SFM_H2O;
                [fluo, sfm] = SFM(data, sfm_wlrange, sfm_H2O);   
                H2O.SFM_SIF = fluo;
                H2O.SFM = sfm;
            elseif strcmp(Name_i(i, end), 'A')   
                sfm_wlrange = windows.SFM.SFM_windowA;
                sfm_O2A = windows.SFM.SFM_O2A;
                [fluo, sfm] = SFM(data, sfm_wlrange, sfm_O2A);    
                O2A.SFM_SIF = fluo;
                O2A.SFM = sfm;
            end
        end
        msg = 'Time delays for SFM: ' + string(toc) + 's';   disp(msg);
    end
    % F_SFM ,by liuxj
    if models.F_SFM == 1 
        disp("Starting F_SFM_"+ string(unique(floor(data.fracDOY))));      tic; 
        [fluo] = F_SFM(data, Flags, Scope);   
        FSFM.FSFM_SIF = fluo';
        msg = 'Time delays for F-SFM: ' + string(toc) + 's';   disp(msg);
    end
    
   % SVD 
    if models.SVD == 1 
        disp("Starting SVD_"+ string(unique(floor(data.fracDOY))));      tic;
        for i = 1:1:4
            Name_i = cell2mat(fieldnames(windows.Others));
            if strcmp(Name_i(i, end), 'H')   
                svd_wlrange = windows.SVD.SVD_windowH;
                svd_Ha = windows.SVD.SVD_Ha;
                rolling = windows.SVD.rolling;                  % 这里在Windows设置
                if ((rolling > 0) && (data.nmeas_ < rolling))
                    % not enough measurements to perform the SVD! NAN
                    SVD = [];
                    SVD.fluo = ones(data.nmeas_, 1) * nan;
                    SVD.refl = ones(data.nmeas_, 1021) * nan;
                else
                    [SVD] = rolling_svd(data, svd_wlrange, svd_Ha, rolling); 
                end 
                Ha.SVD_SIF = SVD.fluo;
                Ha.SVD = SVD;
            elseif strcmp(Name_i(i, end), 'B')   
                svd_wlrange = windows.SVD.SVD_windowB;
                svd_O2B = windows.SVD.SVD_O2B;
                rolling = windows.SVD.rolling;
                if ((rolling > 0) && (data.nmeas_ < rolling))
                    % not enough measurements to perform the SVD! NAN
                    SVD = [];
                    SVD.fluo = ones(data.nmeas_, 1) * nan;
                    SVD.refl = ones(data.nmeas_, 1021) * nan;
                else
                    [SVD] = rolling_svd(data, svd_wlrange, svd_O2B, rolling); 
                end 
                O2B.SVD_SIF = SVD.fluo;
                O2B.SVD = SVD;
            elseif strcmp(Name_i(i, end), 'W')   
                svd_wlrange = windows.SVD.SVD_windowW;
                svd_H2O = windows.SVD.SVD_H2O;
                rolling = windows.SVD.rolling;
                if ((rolling > 0) && (data.nmeas_ < rolling))
                    % not enough measurements to perform the SVD! NAN
                    SVD = [];
                    SVD.fluo = ones(data.nmeas_, 1) * nan;
                    SVD.refl = ones(data.nmeas_, 1021) * nan;
                else
                    [SVD] = rolling_svd(data, svd_wlrange, svd_H2O, rolling); 
                end 
                H2O.SVD_SIF = SVD.fluo;
                H2O.SVD = SVD;
            elseif strcmp(Name_i(i, end), 'A')   
                svd_wlrange = windows.SVD.SVD_windowA;
                svd_O2A = windows.SVD.SVD_O2A;
                rolling = windows.SVD.rolling;
                if ((rolling > 0) && (data.nmeas_ < rolling))
                    % not enough measurements to perform the SVD! NAN
                    SVD = [];
                    SVD.fluo = ones(data.nmeas_, 1) * nan;
                    SVD.refl = ones(data.nmeas_, 1021) * nan;
                else
                    [SVD] = rolling_svd(data, svd_wlrange, svd_O2A, rolling); 
                end 
                O2A.SVD_SIF = SVD.fluo;
                O2A.SVD = SVD;
            end
        end
        msg = 'Time delays for SVD: ' + string(toc) + 's';   disp(msg);
    end
    % DOAS
    if models.DOAS == 1 
        disp("Starting DOAS_"+ string(unique(floor(data.fracDOY))));      tic;
        for i = 1:1:4
            Name_i = cell2mat(fieldnames(windows.Others));
            if strcmp(Name_i(i, end), 'H')   
                doas_wlrange = windows.Others.WindowH;  
                [fluo, ~] = doasFit(data, doas_wlrange);
                Ha.DOAS_SIF = abs(fluo.SIF);
                Ha.DOAS = fluo;
            elseif strcmp(Name_i(i, end), 'B')   
                doas_wlrange = windows.Others.WindowB;  
                [fluo, ~] = doasFit(data, doas_wlrange);
                O2B.DOAS_SIF = abs(fluo.SIF);
                O2B.DOAS = fluo;
            elseif strcmp(Name_i(i, end), 'W')   
                doas_wlrange = windows.Others.WindowW;  
                [fluo, ~] = doasFit(data, doas_wlrange);
                H2O.DOAS_SIF = abs(fluo.SIF);
                H2O.DOAS = fluo;
            elseif strcmp(Name_i(i, end), 'A')   
                doas_wlrange = windows.Others.WindowA;  
                [fluo, ~] = doasFit(data, doas_wlrange);
                O2A.DOAS_SIF = abs(fluo.SIF);
                O2A.DOAS = fluo;
            end
        end
        msg = 'Time delays for DOAS: ' + string(toc) + 's';   disp(msg);
    end
    
    result.Ha = Ha;
    result.O2B = O2B;
    result.H2O = H2O;
    result.O2A = O2A;
    result.FSFM = FSFM';
end

%% -------------------------------------------具体实现------------------------------------ %
% ---sFLD---,Standard Fraunhofer Line Discrimination method (sFLD) from Plascyk and Gabriel 1975, Damm et al. 2011
function [fluo, refl] = sFLD(data, outer, wl_range) 

    nmeas_ = data.nmeas; 
    sky_spec = data.sky';
    veg_spec = data.veg';
    wl = data.wl;
    Eidxo  = wl > outer(1) & wl < outer(2);
    Lidxo  = wl > outer(1) & wl < outer(2); 
    [~, Lidxi, Eidxi] = Cal_Index(wl_range, data);
    
    fluo = zeros(nmeas_, 1);
    refl = zeros(nmeas_, 1);
    for i = 1:nmeas_
        sky_in = sky_spec(i, Eidxi);
        sky_out = mean(sky_spec(i, Eidxo));
        veg_in = veg_spec(i, Lidxi);
        veg_out = mean(veg_spec(i, Lidxo));
        fluo(i) = (sky_out.*veg_in - sky_in.*veg_out)./(sky_out - sky_in); 
        refl(i) = (veg_in - fluo(i)) ./ sky_in;
    end
end

% ---3FLD---, Modified Fraunhofer Line Discrimination method (3FLD) from Maier et al. 2003, Damm et al. 2011
function [fluo, refl] = FLD3(data, left, right, wl_range)  
    nmeas_ = data.nmeas; 
    sky_spec = data.sky';
    veg_spec = data.veg';
    wl = data.wl;

    EidxL  = wl > left(1) & wl < left(2);       % outside / left shoulder
    EidxR = wl > right(1) & wl < right(2);      % right shoulder
    LidxL  = wl > left(1) & wl < left(2);       % outside / left shoulder
    LidxR = wl > right(1) & wl < right(2);      % right shoulder
    [~, Lidxi, Eidxi] = Cal_Index(wl_range, data);

    fluo = zeros(nmeas_, 1);
    refl = zeros(nmeas_, 1);
    for i = 1:nmeas_
        sky_in = sky_spec(i,Eidxi);
        sky_left = mean(sky_spec(i,EidxL));
        sky_right = mean(sky_spec(i,EidxR));
        veg_in = veg_spec(i,Lidxi);
        veg_left = mean(veg_spec(i,LidxL));
        veg_right = mean(veg_spec(i,LidxR));
        wL = (mean(wl(EidxR))-mean(wl(Eidxi)))./(mean(wl(EidxR))-mean(wl(EidxL)));
        wR = (mean(wl(Eidxi))-mean(wl(EidxL)))./(mean(wl(EidxR))-mean(wl(EidxL)));
        fluo(i) = (veg_in - (sky_in./((wL.*sky_left) + (wR.*sky_right))) .* ((wL.*veg_left) + (wR.*veg_right))) ./ (1-(sky_in ./ ((wL.*sky_left) + (wR.*sky_right))));
        refl(i) = (veg_in - fluo(i)) / sky_in;
    end
end

% ---iFLD---,  Improved Fraunhofer Line Discrimination method (iFLD) from Alonso et al. 2008
function [fluo, refl] = iFLD(data, outer, wl_range, Flags)
    nmeas_ = data.nmeas; 
    sky_spec = data.sky';
    veg_spec = data.veg';
    wl = data.wl;
   
    Eidxo  = wl > outer(1) & wl < outer(2);
    Lidxo  = wl > outer(1) & wl < outer(2); 
    [~, Lidxi, Eidxi] = Cal_Index(wl_range, data);
    
    fluo = zeros(nmeas_, 1);
    refl = zeros(nmeas_, 1);
    for i = 1:nmeas_
        % 将吸收线内的波段去掉, 目标——消除吸收线处由荧光填充引起的尖峰
        wl0 = wl;                                     wl0(Flags.Flag == 1) = [];
        sky0 = sky_spec(i, :)';                       sky0(Flags.Flag == 1) = [];
        Ref0 = veg_spec(i, :)' ./ sky_spec(i, :)';    Ref0(Flags.Flag == 1) = [];
        sky_simu(i, :) = interp1(wl0, sky0, wl, 'PCHIP')';
        ref_simu(i, :) = interp1(wl0, Ref0, wl, 'spline')';
        
        sky_in = sky_spec(i, Eidxi);
        sky_out = mean(sky_spec(i, Eidxo));
        veg_in = veg_spec(i, Lidxi);
        veg_out = mean(veg_spec(i, Lidxo));
        Ratio_ref_i = mean(ref_simu(i, Eidxo)) ./ ref_simu(i, Eidxi);
        Ratio_SIF_i = mean(sky_simu(i, Eidxo)) ./ sky_simu(i, Eidxi) * Ratio_ref_i;
        fluo(i)  = (Ratio_ref_i .* sky_out .* veg_in - sky_in * veg_out) ...
                    ./ (Ratio_ref_i * sky_out - Ratio_SIF_i * sky_in);
        refl(i) = (veg_in - fluo(i))./ sky_in;
    end
    

end

% ---pFLD---, by liuxj
function [fluo, refl] = pFLD(data, outer, wl_range, Flags, Scope)
    nmeas_ = data.nmeas; 
    sky_spec = data.sky';
    veg_spec = data.veg';
    wl = data.wl;
   
    Eidxo  = wl > outer(1) & wl < outer(2);
    Lidxo  = wl > outer(1) & wl < outer(2); 
    [~, Lidxi, Eidxi] = Cal_Index(wl_range, data);
    
    fluo = zeros(nmeas_, 1);
    refl = zeros(nmeas_, 1);
    for i = 1:nmeas_
        % 将吸收线内的波段去掉, 目标——消除吸收线处由荧光填充引起的尖峰
        wl0 = wl;                                     wl0(Flags.Flag == 1) = [];
        sky0 = sky_spec(i, :)';                       sky0(Flags.Flag == 1) = [];
        Ref = veg_spec(i, :)' ./ sky_spec(i, :)';
        sky_simu(i, :) = interp1(wl0, sky0, wl, 'PCHIP')';

        [A_p, ~, ~] = lscov(Scope.PC8_Ref, Ref, 1 - Flags.Flag);
        ref_simu(i, :) = Scope.PC8_Ref * A_p;
        
        sky_in = sky_spec(i, Eidxi);
        sky_out = mean(sky_spec(i, Eidxo));
        veg_in = veg_spec(i, Lidxi);
        veg_out = mean(veg_spec(i, Lidxo));
        Ratio_ref_i = mean(ref_simu(i, Eidxo)) ./ ref_simu(i, Eidxi);
        Ratio_SIF_i = mean(sky_simu(i, Eidxo)) ./ sky_simu(i, Eidxi) * Ratio_ref_i;
        fluo(i)  = (Ratio_ref_i .* sky_out .* veg_in - sky_in * veg_out) ...
                    ./ (Ratio_ref_i * sky_out - Ratio_SIF_i * sky_in);
        refl(i) = (veg_in - fluo(i))./ sky_in;
    end
end

% ---SFM---,Spectral Fitting Method from Meroni et al. 2009
function [fluo, sfm] = SFM(data, sfm_wlrange, sfm_Band)
    nmeas_ = data.nmeas; 
    sky_spec = data.sky';
    veg_spec = data.veg';
    wl = data.wl;
    
    idxr = find(wl > min(sfm_wlrange) & wl < max(sfm_wlrange));
    Range = find(wl > sfm_Band(1) & wl < sfm_Band(2));
    Lidxi = Range(1);
    idx = wl(idxr) == wl(Lidxi);
    
    sfm = [];
    sfm.refl = zeros(nmeas_, length(idxr));
    sfm.reduced_chisq = zeros(nmeas_, 1);
    sfm.rmse = zeros(nmeas_, 1);
    sfm.rrmse = zeros(nmeas_, 1);
    sfm.spec_meas = zeros(nmeas_, length(idxr))*NaN;
    sfm.spec_mod = zeros(nmeas_, length(idxr))*NaN;
    for i = 1:nmeas_
        dlambda = wl(idxr) - wl(idxr(1));
        p1 = dlambda.^2 .* sky_spec(i, idxr)';
        p2 = dlambda .* sky_spec(i, idxr)';
        p3 = sky_spec(i, idxr)';
        p4 = dlambda.^2;
        p5 = dlambda;
        p6 = ones(size(idxr));
        K = [p1, p2, p3, p4, p5, p6];
        b = K\veg_spec(i, idxr)';
        y_hat = K*b;
        y = veg_spec(i, idxr)';
        rmse(i) = sqrt(mean((y_hat - y).^2));
        rrmse(i) = rmse(i) ./ mean(y);
        fluo_ = b(4).*dlambda.^2 + b(5).*dlambda + b(6);
        refl(i,:) = (veg_spec(i,idxr) - fluo_') ./ sky_spec(i,idxr);
        fluo(i) = fluo_(idx);
        
        nparams = size(K,2);
        sv_end = b(end);
        [~, sfm.reduced_chisq(i)] = reduced_chisq_screen(y', y_hat', length(idxr), nparams, sv_end, 1);
        [sfm.rmse(i), sfm.rrmse(i), ~] = error_calc(y', y_hat');
    
        sfm.spec_meas(i,:) = y';
        sfm.spec_mod(i,:) = y_hat';
    end
    [fluo] = quality_filter(fluo, sfm.reduced_chisq);
    fluo = fluo';
    sfm.refl = refl;
    sfm.wavewindow = wl(idxr);
end

% ---F_SFM---, by liuxj
function [fluo] = F_SFM(data, Flags, Scope)  
    nmeas_ = data.nmeas; 
    sky_spec = data.sky';
    veg_spec = data.veg';
    
    Reflected = @(a,x)a(1).*x(:,1) + a(2).*x(:,2) + a(3).*x(:,3) + a(4).*x(:,4) + a(5).*x(:,5) + a(6).*x(:,6) + a(7).*x(:,7) + a(8).*x(:,8);
    Fluor_PCsimu = @(a,x)a(9).*x(:,9) + a(10).*x(:,10) + a(11).*x(:,11) + a(12).*x(:,12) + a(13).*x(:,13);
    TotalFit = @(a,x)Reflected(a,x) + Fluor_PCsimu(a,x);
    F_result_PCA = 0;
    
    for i = 1:nmeas_
        for m = 1:1:3
            Ref_ratio_f = (veg_spec(i, :)' - F_result_PCA) ./ sky_spec(i, :)';
            
            [A_f, ~, ~] = lscov(Scope.PC8_Ref, Ref_ratio_f, 1 - Flags.Flag);
            Ref_simu_SFM = Scope.PC8_Ref * A_f;
            
            RadRef_H = sky_spec(i, :)' .* Ref_simu_SFM .* Flags.Flag_H;
            RadRef_B = sky_spec(i, :)' .* Ref_simu_SFM .* Flags.Flag_B;
            RadRef_W = sky_spec(i, :)' .* Ref_simu_SFM .* Flags.Flag_W;
            RadRef_A = sky_spec(i, :)' .* Ref_simu_SFM .* Flags.Flag_A;
            Rad_H = sky_spec(i, :)' .* Flags.Flag_H;
            Rad_B = sky_spec(i, :)' .* Flags.Flag_B;
            Rad_W = sky_spec(i, :)' .* Flags.Flag_W;
            Rad_A = sky_spec(i, :)' .* Flags.Flag_A;
            
            QE65_veg_used = veg_spec(i, :)';
            QE65_veg_used(Flags.Flag <= 0) = [];
            XX = [RadRef_H, Rad_H, RadRef_B, Rad_B, RadRef_W, Rad_W, RadRef_A, Rad_A, Scope.PC5_Fluor];
            X = XX;
            X(Flags.Flag <= 0,:) = [];
            
            a = [1 0 1 0 1 0 1 0 0 0 0 0 0];
            lb = [ 0 -1  0 -1  0 -1  0 -0.5     0  -0.05 -0.01 -0.01 -0.05];
            ub = [ 2  1  2  0  2  1  2  0.5  0.05   0.05  0.01  0.01  0.05];
            [A0, ~] = lsqcurvefit(TotalFit,a,X,QE65_veg_used,lb,ub);
            F_result_PCA = Fluor_PCsimu(A0,XX);
        end
        fluo(:, i) = F_result_PCA;
    end
    
end

% ---SVD---,Singular vector decomposition method, inspired by code from Ari Kornfeld and Christian Frankenberg
function [svd_] = rolling_svd(data, svd_wlrange, svd_Band, rolling)   
% Also includes automated optimization, see below.
% *** Step 1: define default parameters ***
    sky_spec = data.sky';
    veg_spec = data.veg';
    wl = data.wl;
    hf = data.hf;
    
    Range = find(wl > svd_Band(1) & wl < svd_Band(2));
    Lidxi = Range(1);
    pixels = wl > min(svd_wlrange) & wl < max(svd_wlrange);
    % prescribed SIF shape (Guanter et al. 2013 RSE)
    Fluo = hf(Lidxi);
    hf_ = hf(pixels);
    
    spectra = sky_spec(:, pixels);
    spec_sky = spectra;
    
    np1_combos = linspace(3,6,4);  
    np2_combos = linspace(3,6,4);
    sv_combos = combvec(1, np1_combos, np2_combos);
    
    np1_ = 0;
    np2_ = 0;
    numSV_ = 0;
    minBIC = 1000000000;  % start with arbitrarily large number...
    svd_ = [];
    svd_.numSV = [];
    svd_.EV = [];
    svd_.rmse = [];
    svd_.rrmse = [];
    svd_.rmsp = [];
    
    norm_ = 1;            % add by liuxj
    
    if rolling > 1        % will use middle of rolling window as spec veg
        [n, ~] = size(spectra);
        spec_meas = zeros(n-rolling+1, length(pixels));
        spec_mod = zeros(n-rolling+1, length(pixels));
        sv_end = zeros(n-rolling+1, 1);
        for j=1:length(sv_combos)
            for k = 1:(n-rolling+1)
                spec_sky_window = spec_sky(k:(k+rolling-1),:);
                [K_frag, ~, ~, numSV, EV, all_V, norm_S, nparams] = train_svd_(spec_sky_window, wl, pixels, hf_, sv_combos(2,j), sv_combos(3,j), norm_);
                [~, spec_mod(k,:), spec_meas(k,:), sv_frag] = solve_svd_(K_frag, veg_spec(k+ceil(rolling/2)-1, :), pixels);
                sv_end(k) = sv_frag(end);
            end
            [minBIC, svd_] = optimize_svd(spec_mod, spec_meas, pixels, Fluo, minBIC, np1_, np2_, numSV_, sv_combos(2,j), sv_combos(3,j), numSV, rolling, EV, norm_S, svd_, nparams, sv_end, sky_spec(1:(n-rolling+1),:), veg_spec(1:(n-rolling+1),:), hf, wl, all_V);
        end
    elseif rolling < 0                           % use x time points as training spectra
        n = abs(rolling);
        midx = (1:n)*floor(size(spectra,1)/n);   % 25 scattered throughout data
        for j=1:length(sv_combos)
            [K, ~, ~, numSV, EV, all_V, norm_S, nparams] = train_svd_(spec_sky(midx,:), wl, pixels, hf_, sv_combos(2,j), sv_combos(3,j), norm_);
            [~, spec_mod, spec_meas, sv] = solve_svd_(K, veg_spec, pixels);
            sv_end = sv(:,end);
            [minBIC, svd_] = optimize_svd(spec_mod, spec_meas, pixels, Fluo, minBIC, np1_, np2_, numSV_, sv_combos(2,j), sv_combos(3,j), numSV, rolling, EV, norm_S, svd_, nparams, sv_end, sky_spec, veg_spec, hf, wl, all_V);
        end
    else
        for j = 1:length(sv_combos)
            [K, ~, ~, numSV, EV, all_V, norm_S, nparams] = train_svd_(spec_sky, wl, pixels, hf_, sv_combos(2,j), sv_combos(3,j), norm_);
            [~, spec_mod, spec_meas, sv] = solve_svd_(K, veg_spec, pixels);
            sv_end = sv(:,end);
            [minBIC, svd_] = optimize_svd(spec_mod, spec_meas, pixels, Fluo, minBIC, np1_, np2_, sv_combos(2,j), sv_combos(3,j), numSV, rolling, EV, norm_S, svd_, nparams, sv_end, sky_spec, veg_spec, hf, wl, all_V);
        end 

    end
end
function [K, Ksd, FLidx, numSV, EV, all_V, norm_S, nparams] = train_svd_(training_spec, wl, pixels, hf, np1, np2, norm_)
    % *** Step 2: train SVD ***
    [~,S,V] = svd(training_spec, 'econ');
    numSV = 20;                              % ***---注意：这里限制,很重要
    norm_S = diag(S)/sum(diag(S))*100;
    
    % Return at least the first 2 Eigenvectors passing the eigenvalue threshold :
    evToreturn = 1:min( max(2, numSV), size(V, 2)); % at least 2
    EV = -V(:,evToreturn);  % Singular vectors always come out -ve, so reverse them here. 
    all_V = -V;
    
    % Compute Matrices for the SV convolved with polynomial terms in lambda space which is needed to avoid spurious correlation of the power terms
    poly1 = zeros(size(EV,1),np1);
    poly2 = zeros(size(EV,1),np2);

    % 1st SV
    idx = 1;
    for i=np1:-1:0 %0
        poly1(:,idx) =  EV(:, 1).*(wl(pixels) - mean(wl(pixels))).^i;
        idx = idx + 1;
    end

    if size(EV,2) > 1
        % 2nd SV
        idx = 1;
        for i=np2:-1:0 %0
            poly2(:,idx) =  EV(:, 2).*(wl(pixels) - mean(wl(pixels))).^i;
            idx = idx + 1;
        end

        % remaining SV
        EV3plus = EV(: , 3:numSV);
        K = [poly1, poly2, EV3plus, hf];
    else
        K = [poly1, hf];
    end
    nparams = size(K,2);
    
    FLidx = (size(K, 2) - size(hf, 2)+1):size(K,2);
    Ksd = ones(1, size(K, 2));
    if norm_ == 1
        %Ari's version: normalize results so SD = 1 but do not center
        Ksd = std(K); %to disable, use: ones(1, size(R.K, 2));
%         Kmean = zeros(1,size(K, 2)); %mean(K); 
%         K = MxV(@minus, K, Kmean);
%         K = MxV(@rdivide, K, Ksd);
    end
end
function [result, spec_mod, spec_meas, sv] = solve_svd_(K, spec_veg, pixels)
% ** Step 3: solve SVD *********************************************
    spec_meas = spec_veg(:,pixels);
    result = K\spec_meas';
    spec_mod = (K*result)'; %modeled spectrum
    sv = result';
end
function [minBIC, svd_] = optimize_svd(spec_mod, spec_meas, pixels, Fluo, minBIC, ~, ~, np1, np2, numSV, rolling, EV, norm_S, svd_, nparams, sv_end, spec_sky, spec_veg, hf, wavelengths, all_V)
% %% Step 4: optimize SVD using BIC ***************************
    norm_spec_mod = zeros(size(spec_mod));
    norm_spec_meas = zeros(size(spec_meas));
    
    %rolling_ = "see prev";
    
    for i=1:size(spec_mod,1)
        norm_spec_mod(i,:) = (spec_mod(i,:) - min(spec_mod(i,:))) / (max(spec_mod(i,:) - min(spec_mod(i,:))));
        norm_spec_meas(i,:) = (spec_meas(i,:) - min(spec_meas(i,:))) / (max(spec_meas(i,:) - min(spec_meas(i,:))));
    end
    
    BIC = zeros(1,size(spec_mod,1));
    AIC = zeros(1,size(spec_mod,1));
    for i=1:size(spec_mod,1)
        resids = norm_spec_meas(i,:)-norm_spec_mod(i,:);
        mu = nanmean(resids);
        sigma = nanstd(resids); 
        L = (1/(sqrt(2*pi)*sigma))^size(resids,2) * exp(-sum(((resids-mu).^2)/2*sigma.^2));
        LL = -(size(resids,2))*log(sqrt(2*pi))-(size(resids,2))*log(sigma) - sum(((resids-mu).^2)/(2*sigma.^2));
        BIC(i) = nparams*log(length(pixels)) - (2.*(LL));
        AIC(i) = 2*nparams - 2.*(LL);
    end
    SIF_range = zeros(length(sv_end), length(hf));
    if nanmean(BIC) < minBIC
        
        [SIF, reduced_chisq] = reduced_chisq_screen(spec_meas, spec_mod, size(spec_meas,2), nparams, sv_end, Fluo);
        [rmse, rrmse, rmsp] = error_calc(spec_meas, spec_mod);
        [SIF2] = quality_filter(SIF, reduced_chisq);
        minBIC = nanmean(BIC);
        np1_ = np1;
        np2_ = np2;
        numSV_ = numSV;
%         rolling_ = rolling;
        
        SIF2 = [zeros(ceil(rolling/2)-1,1) * NaN; SIF2; zeros(ceil(rolling/2)-1,1) * NaN]; %kludge to insert NaNs if using rollingsvd
        for i = 1:length(sv_end)
            SIF_range(i,:) = sv_end(i) .* hf;
        end
        svd_.refl = (spec_veg-SIF_range) ./ spec_sky; 
        svd_.refl = [zeros(ceil(rolling/2)-1,size(svd_.refl,2))*NaN; svd_.refl; zeros(ceil(rolling/2)-1,size(svd_.refl,2))*NaN];
        svd_.spec_meas = spec_meas;
        svd_.spec_mod = spec_mod;
        svd_.wavewindow = wavelengths(pixels);
        svd_.fluo = SIF2;
        svd_.numSV = [svd_.numSV, numSV_];
        svd_.EV = EV;
        svd_.all_V = all_V;
        svd_.norm_S = norm_S;
        svd_.reduced_chisq = reduced_chisq;
        svd_.rmse = [svd_.rmse, rmse];
        svd_.rrmse = [svd_.rrmse, rrmse];
        svd_.rmsp = [svd_.rmsp, rmsp];
        svd_.np1 = np1_;
        svd_.np2 = np2_;
    end
    %msg = "BIC: "+string(minBIC)+", rolling = "+rolling_+"; numSV = ["+string(min(svd_.numSV))+"-"+string(max(svd_.numSV))+"]; np1 = "+string(np1_)+"; np2 = "+string(np2_);
    %disp(msg);
end

% ---DOAS---,Differential Optical Absorption Spectroscopy method (Grossmann et al. 2018) modified from Christian Frankenberg's scripts
 function [FR,xx] = doasFit(data, doas_wlrange) 
    nmeas_ = data.nmeas; 
    sky_spec = data.sky';
    veg_spec = data.veg';
    wl = data.wl;
    hf = data.hf;
    
    [~, Lidxi, ~] = Cal_Index(doas_wlrange, data);
    idxr = find(wl > min(doas_wlrange) & wl < max(doas_wlrange));
    
	xx = ((wl(idxr))-mean(wl(idxr)))/100; 

    hf_meanctr = (hf/mean(hf));
    h = hf_meanctr(idxr);
    Fluo_in = hf_meanctr(Lidxi);

    xl = xx/max(abs(xx)); %normalize wavelengths into -1 to 1 range
    
	nn = length(idxr);
    le = nmeas_;
    FR.SIF = zeros(le,1)*NaN;
    FR.reduced_chisq = zeros(le,1)*NaN;
    FR.rmse = zeros(le,1)*NaN;
    FR.rrmse = zeros(le,1)*NaN;
    FR.SIF_1sigma = zeros(le,1)*NaN;
    FR.res = zeros(le,nn)*NaN;
    FR.spec_meas = zeros(le,nn)*NaN;
    FR.spec_mod = zeros(le,nn)*NaN;

    K = [];
	legendre_num = 6; 
    
    for j = 0:legendre_num 
        K = [K legendreP(j,xl)];
    end

    FR.refl = zeros(le, length(wl))*NaN;

    for i=1:nmeas_

        diffSVD = [K h./veg_spec(i,idxr)']; %Philipp Kohler's Legendre polynomial solution

        y = log(veg_spec(i,idxr))-log(sky_spec(i,idxr));
        x = diffSVD\y';

        mod = diffSVD*x;
        resid = mod-y';

		Ssquared = (mod-y').^2; 
		FR.res(i,:) = resid; 
		FR.spec_meas(i,:) = y;
		FR.spec_mod(i,:) = mod;
		FR.wavewindow = wl(idxr);
		S = sum(Ssquared)/(length(y)-length(x))*inv(transpose(diffSVD)*diffSVD);
		FR.SIF_1sigma(i) = sqrt(S(end,end));
		FR.idxr = idxr;
		nparams = size(diffSVD,2);
		sv_end = x(end);
		[FR.SIF(i), FR.reduced_chisq(i)] = reduced_chisq_screen(y, mod', length(idxr), nparams, sv_end, Fluo_in);
		SIF_range = sv_end .* hf_meanctr;
		FR.refl(i,:) = (veg_spec(i,:)-SIF_range') ./ sky_spec(i,:);
		[~, FR.rrmse(i), ~] = error_calc(y, mod');
    end
    
    [FR.SIF] = quality_filter(FR.SIF, FR.reduced_chisq);
end


% ---Support functions---,
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
% --函数备注：
function [SIF, reduced_chisq] = reduced_chisq_screen(spec_meas, spec_mod, npx, nparams, sv_end, Fluo_in)
% calculate reduced chi squared, 指简化的卡方检验,RSS/dof 残差均方,标准化残差均方根（standardized root mean square residual，SRMR）
    degf = npx - nparams; 
    reduced_chisq = zeros(size(spec_meas,1),1)*NaN;
    SIF = sv_end * Fluo_in;
    for i = 1:size(spec_meas,1)
        var_resid = var(spec_meas(i,:) - spec_mod(i,:), 0, 2);           % 沿第二个维度计算的方差。
        chisq = nansum((spec_meas(i,:)- spec_mod(i,:)).^2 ./ var_resid); 
        reduced_chisq(i) = chisq/degf;
    end
end
% --函数备注：
function [SIF] = quality_filter(SIF, reduced_chisq)
% Check for goodness of fit using reduced chi squared (Guanter et al. 2012, Sun et al. 2018)
    chisq_max = 2;
    chisq_min = 0.5;
    for i=1:size(SIF,1)
        if (reduced_chisq(i) > chisq_max) || (reduced_chisq(i) < chisq_min) 
            SIF(i) = NaN;
        end
    end
end
% --函数备注：
function [rmse, rrmse, rmsp] = error_calc(spec_meas, spec_mod)

    rmse = sqrt(nanmean((spec_mod - spec_meas).^2, 2));
    rrmse = (rmse ./ abs(nanmean(spec_meas,2)) .* 100);
    rmsp = rmse ./ spec_meas .* 100;
end
