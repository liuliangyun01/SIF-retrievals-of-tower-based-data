
function SZA = Cal_SZA(time_fixed, year, month, day, lon, lat, TimeZone )

    for t = 1:length(time_fixed)
        t_i = time_fixed(t);
        hour = floor(t_i*24);
        min = floor((t_i*24-hour)*60);
        sec = floor(((t_i*24-hour)*60-min)*60);
        % 年积日的计算
        JD1 = fix(365.25*(year-1))+fix(30.6001*(1+13))+1+hour/24+1720981.5;
        if month<=2
            JD2 = fix(365.25*(year-1))+fix(30.6001*(month+13))+day+hour/24+1720981.5;
        else
            JD2 = fix(365.25*year)+fix(30.6001*(month+1))+day+hour/24+1720981.5;
        end

        DOY = JD2-JD1+1;
        N0 = 79.6764 + 0.2422*(year-1985) - floor((year-1985)/4.0);
        sitar = 2*pi*(DOY-N0)/365.2422;
        ED = 0.3723 + 23.2567*sin(sitar) + 0.1149*sin(2*sitar) - 0.1712*sin(3*sitar)- 0.758*cos(sitar) + 0.3656*cos(2*sitar) + 0.0201*cos(3*sitar);
        ED = ED*pi/180;                   % //ED本身有符号，无需判断正负。
        if (TimeZone == -13)
            dLon = lon - (floor((lon*10-75)/150)+1)*15.0;
        else
            dLon = lon - TimeZone*15.0;   % //地球上某一点与其所在时区中心的经度差
        end
        Et = 0.0028 - 1.9857*sin(sitar) + 9.9059*sin(2*sitar) - 7.0924*cos(sitar)- 0.6882*cos(2*sitar); % //视差
        gtdt = hour + min/60.0 + sec/3600.0 + dLon/15;  % //地方时
        gtdt = gtdt + Et/60.0;
        dTimeAngle = 15.0*(gtdt-12);
        dTimeAngle = dTimeAngle*pi/180;
        latitudeArc = lat*pi/180;
        HeightAngleArc = asin(sin(latitudeArc)*sin(ED)+cos(latitudeArc)*cos(ED)*cos(dTimeAngle));       % 计算太阳高度角
        SZA(t) = 90.0 - HeightAngleArc*180/pi;
    end
end
     
     