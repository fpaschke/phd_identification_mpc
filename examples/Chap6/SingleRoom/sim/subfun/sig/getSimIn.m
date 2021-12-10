function In = getSimIn(P)
% Return inputs for simulation. Will generate inputsignals for one week
% more since MPC needs predictions.
%
% In = getSimIn(P)
% In [struct]:  struct with timeseries data necessry for simulation
% P [struct]:   parameter struct necessary for generation of inputs. See getSimPar()

In.Time = (datenum(P.Tstart) - datenum([P.Tstart(1:4) '-1-1']))*86400+(0:P.T_s*60:(P.N_w+1)*7*60*24*60)'; % Time in s
if ~P.Sig.N_pers.Periodic
    In.N_pers = timeseries(buildRndOccSig(In.Time,P.Sig.N_pers.NumPeop,P.Sig.N_pers.Schedule),In.Time); % Occupancy
else
    In.N_pers = timeseries(buildPerOccSig(P.N_w+1,P.T_s,P.Sig.N_pers.NumPeop,P.Sig.N_pers.StdPeop,P.Sig.N_pers.Schedule,P.Sig.N_pers.RiseTime),In.Time); % Occupancy
end
% Environment
In.WeatherData = importdata(P.WeatherFile);
In.T_out = timeseries(interp1([In.WeatherData.data(:,1); In.WeatherData.data(2:end,1)+In.WeatherData.data(end,1)],[In.WeatherData.data(:,2); In.WeatherData.data(2:end,2)],In.Time),In.Time,'Name','Outside Temp.'); % Outside Temp
In.P_sun = timeseries(interp1([In.WeatherData.data(:,1); In.WeatherData.data(2:end,1)+In.WeatherData.data(end,1)],[In.WeatherData.data(:,9); In.WeatherData.data(2:end,9)],In.Time),In.Time,'Name','Global hor rad.'); % Global hor rad. [Wh/m2]
%In.P_sun_dirnor = timeseries(interp1(In.WeatherData.data(:,1),In.WeatherData.data(:,10),In.Time),In.Time,'Name','Globalstrahlung')/1e3; % Global hor rad. [kWh/m2]
%In.P_sun_diff = timeseries(interp1(In.WeatherData.data(:,1),In.WeatherData.data(:,11),In.Time),In.Time,'Name','Diffuse Strahlung')/1e3; % Global hor rad. [kWh/m2]
[El, Az] = util.calcSunPosition(In.Time/86400+datenum([P.Tstart(1:4) '-1-1']),P.GeoPos(2),P.GeoPos(1));
In.Azimuth = timeseries(Az,In.Time,'Name','Azimutwinkel'); % Outside Temp
In.Elevation = timeseries(El,In.Time,'Name','Sonnenhöhe'); % Outside Temp

% Build Upper and lower Temp. Setpoint signals
[Tl,Tu] = buildTsetSig(In.Time,P.Sig.T_ref,In.T_out.Data);                     
In.T_ref_lower = timeseries(Tl,In.Time,'Name','Lower temp. bound'); 
In.T_ref_upper = timeseries(Tu,In.Time,'Name','Upper temp. bound'); 

% Build Identification Signals
In.H_flo = timeseries(buildRndSig(P.N_w+1,P.T_s,P.Sig.H_flo,0),In.Time,'Name','Heizventilposition'); 
In.H_win = timeseries(buildRndSig(P.N_w+1,P.T_s,P.Sig.H_win,0),In.Time,'Name','Verschattungsposition'); 

% Set Time Properties
In = setTimeProp(In,[P.Tstart(1:4) '-1-1'],P.DateFormat);
end