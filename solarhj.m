function [maxwattstable, minwattstable] = solarhj;
%%% SOLARHJ is a function to compute the minimum and maximum solar radiation possible in W/m2 from the stations 
% at the HJ Andrews. 

% the outputs are 2 "tables" containing the array of values representing the max and minimum W/m2 that are available 
% "instantaneously" in the order of 'PRIMET','CS2MET','CENMET','VANMET','H15MET', and 'UPLMET'.

% this function can be used within another QAQC script; call 
% [max_table, min_table] = solarhj;
% this will enable you to have the look up tables in another file for reference.

% Created by Fox Peterson
% o7-31-2014
% SUN_POSITION is a BIF from MATLAB
% the algorithm used to compute irradiance is from the help manual for ARC GIS 10.2, solar insolation tool
% however, the Transmissivity Parameter comes from the construction for Air Mass presented in the 
% online lectures from deAnza College at www.deanza.edu/faculty/hamidiridha. The algorithm has been
% altered by Fox
% the cloud cover constructs come from the national solar radiation db handbook 1961-1990, 
% which is available online at rredc.nrel.gov
% the calculation of incidence angle comes from notes taken by Fox Peterson while working with Chris and Sean 
% Carigg, which I have confirmed against the python code in the google earth engine, see code.google.com.

    SCONST       = 1367;
    TIMEINT      = 1;     % hours (1 hour intervals here...)
    cloudcover_h = 0;
    cloudcover_m = 0;
    cloudcover_l = 0;
    utc_offset   = -8;
    M_theta      = 1;    % Relative optical length, representative of "air mass"; should be between 0 and 2.5.

    % a small table containing the values of the stations to put into sun pos
    stations = {{'PRIMET',{[44.21189300, -122.25594100], 436, 4.80, 158}}; 
                {'CS2MET',{[44.21473500, -122.24862900], 482, 10.26, 19}};
                {'CENMET',{[44.24348200, -122.14109300], 1028, 12.17, 224}};
                {'VANMET',{[44.27161800,-122.14936800], 1268,5.48, 183}};
                {'H15MET',{[44.2645069, -122.173778247], 909, 29.60, 215}};
                {'UPLMET',{[44.20709700, -122.11976300], 1298, 8.31, 72}}
                };

    % a time series ranging from 2010-01-01 00:00:00 to 2010-12-31 23:59:59 by 1 hour breaks           
    timeseries1 = [734139:0.0416666667442769:734503.958333333];

    Qstar_stations = zeros(length(timeseries1),6);
    Qstar = zeros(length(timeseries1),1);

    for i = 1:len(stations)
        stationname        = stations{i}{1} 
        location.latitude  = stations{i}{2}{1}(1);
        location.longitude = stations{i}{2}{1}(2);
        location.altitude  = stations{i}{2}{2};
        stationslope       = stations{i}{2}{3};
        stationaspect      = stations{i}{2}{4};
        
        for j = 1:length(timeseries1)
                datetime    =   datevec(timeseries1(j));
                time.year   =   datetime(1,1);       % Valid for [-2000, 6000]
                time.month  =   datetime(1,2);       % month [1-12]
                time.day    =   datetime(1,3);       % calendar day [1-31]
                time.hour   =   datetime(1,4);       % local hour [0-23]
                time.min    =   datetime(1,5);       % minute [0-59]
                time.sec    =   datetime(1,6);       % second [0-59]
                time.UTC    =   utc_offset; 
                sun = sun_position(time,location);   % sun.zenith, sun.azimuth
                solarazimuth = sun.azimuth;
                solarzenith  = sun.zenith;

                % compute the parameters needed for solar insolation

                %%% function here assumes we have a perfectly clear atmosphere but that earth is spherical
                if solarzenith > 0 && solarzenith < 180;
                    Transmissivity = 1/(4*pi()^2)*sind(solarzenith)*(1-0.7*cloudcover_h)*(1-0.4*cloudcover_m)*(1-0.4*cloudcover_l);
                else 
                    Transmissivity == 0.1;
                end

                AngleIncidence = cosd(solarazimuth - stationaspect)*sind(stationslope)*sind(solarzenith) - cosd(solarzenith)*cosd(stationslope);
                %HillEFX = @(solarazimuth, solarzenith, stationaspect, stationslope) cos(rad2deg(solarazimuth) - rad2deg(stationaspect)) * sin(rad2deg(stationslope))*sin(rad2deg(solarzenith)) + cos(rad2deg(solarzenith))*cos(rad2deg(stationslope));

                Qstar(j,1) = SCONST * Transmissivity *24* cosd(AngleIncidence);
        end

        %clearvars stationname location stationslope station aspect

        Qstar_stations(:,i) = Qstar;
        %Qstar = SWup + SWdown + LWup + LWdown;
    end

    maxwattstable = Qstar_stations;

    filename = 'Maximum_solar_downwelling_Wm2.csv';
    fid = fopen(filename,'w');
    fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n','DATETIME','PRIMET','CS2MET','CENMET','VANMET','H15MET','UPLMET');
    
    [rows,cols] = size(Qstar_stations)
    for i = 1:rows
    fprintf(fid,'%s,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f',datestr(timeseries1(i)),maxwattstable(i,1:cols-1));
    fprintf(fid,'%7.3f\n',Qstar_stations(i,cols));
    end
    fclose(fid);

    % assume cloudcover_h = 1, cloudcover_l = 1, and cloudcover_m = 1-- fully clouds
    % now assume that this means the diffuse removal is then (0.3)*(0.6)*(0.6) = 0.1080
    minwattstable = Qstar_stations*0.1080;

    filename = 'Miniumum_solar_downwelling_Wm2.csv';
    fid = fopen(filename,'w');
    fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n','DATETIME','PRIMET','CS2MET','CENMET','VANMET','H15MET','UPLMET');
    
    [rows,cols] = size(Qstar_stations)
    for i = 1:rows
    fprintf(fid,'%s,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f',datestr(timeseries1(i)), minwattstable(i,1:cols-1));
    fprintf(fid,'%7.3f\n',minwattstable(i,cols));
    end
    fclose(fid);

end
