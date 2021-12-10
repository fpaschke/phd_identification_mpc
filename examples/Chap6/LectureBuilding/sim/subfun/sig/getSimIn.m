function In = getSimIn(P)
% Return inputs for simulation. Will generate inputsignals for one week
% more since MPC needs predictions.
%
% In = getSimIn(P)
% In [struct]:  struct with timeseries data necessry for simulation
% P [struct]:   parameter struct necessary for generation of inputs. See getSimPar()

In.Time = (datenum(P.Tstart) - datenum([P.Tstart(1:4) '-1-1']))*86400+(0:P.T_s*60:(P.N_w+1)*7*60*24*60)'; % Time in s
In.N_pers = timeseries(cell2mat(arrayfun(@(nh) buildPerOccSig(P.N_w+1,P.T_s,squeeze(P.Sig.N_pers.Mean(nh,:,:)),P.Sig.N_pers.Std,...
                        P.Sig.N_pers.Schedule,P.Sig.N_pers.T_r,nh),1:size(P.Sig.N_pers.Mean,1),'UniformOutput',0)),In.Time); % Occupancy

% Environment
In.WeatherData = importdata(P.WeatherFile);
In.T_out = timeseries(interp1([In.WeatherData.data(:,1); In.WeatherData.data(2:end,1)+In.WeatherData.data(end,1)],[In.WeatherData.data(:,2); In.WeatherData.data(2:end,2)],In.Time),In.Time,'Name','Outside Temp.'); % Outside Temp
In.P_sun = timeseries(interp1([In.WeatherData.data(:,1); In.WeatherData.data(2:end,1)+In.WeatherData.data(end,1)],[In.WeatherData.data(:,9); In.WeatherData.data(2:end,9)],In.Time),In.Time,'Name','Global hor rad.'); % Global hor rad. [Wh/m2]
[El, Az] = util.calcSunPosition(In.Time/86400+datenum([P.Tstart(1:4) '-1-1']),P.GeoPos(2),P.GeoPos(1));
In.Azimuth = timeseries(Az,In.Time,'Name','Sun Azimuth'); % Outside Temp
In.Elevation = timeseries(El,In.Time,'Name','T_out'); % Outside Temp
              
% Build Upper and lower Temp. Setpoint signals
[Tl, Tu] = buildTsetSig(In.Time,P.Sig.T_ref,In.T_out.Data);                     
In.T_ref_lower = timeseries(Tl,In.Time,'Name','Temp. bounds'); 
In.T_ref_upper = timeseries(Tu,In.Time,'Name','Temp. bounds');  

% Build Identification Signals
for nz = 1:length(P.Sig.m_sup)
    m_sup(:,nz) = buildRndSig(P.N_w,P.T_s,P.Sig.m_sup(nz),nz);
end
N = length(m_sup);
In.m_sup = timeseries(m_sup,In.Time(1:N),'Name','Zuluftvolumenstrom'); 
In.T_sup = timeseries(buildRndSig(P.N_w,P.T_s,P.Sig.T_sup,nz+1),In.Time(1:N),'Name','Zulufttemperatur'); 

In = setTimeProp(In,[P.Tstart(1:4) '-1-1'],P.DateFormat);
end