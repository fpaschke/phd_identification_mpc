function M = estWeatherModel(Data,P)
% Setups the Model structure, estimates models from data
% and displays the identification result. At the end user will be prompted
% to select model orders. See also: idModels.NsfPolyModel.estModel for more
% information on parameters.
%
% M [NsfPolyModel]:         Estimated model
% Data [struct or object]:  Contains timeseries data necessary for identification (N #Samples):
%                           All Signals are assumed to be uniformly sampled.
%           Data.T_out [1 x N]: Outside Air Temperature [C] (input)
%           Data.P_sun [1 x N]: Global Radiation on horizontal plane [kW/m2] (input)
% P [struct]: Structure containing Model properties validation and identification options.
%           P.fc [double]:          Cutoff Frequency of 1st order highpass butterworth prefilter in [1/h]. 
%                                   1/(6*2*pi) corresponds to 1st Order System with Timeconstant of 6h
%           P.Orders [int vector]:  Orders of Polynomials A,B,... for which identification will be performed
%           P.Seasonality (int or []): Seasonality of NoiseModel in number of samples
%           P.SeasonalityOrder [logical]: Seasonality will be Increased by Order if 'ByOrder'
%           P.HpVal [int]:          Prediction horizon in #steps for validation (plot)
%           P.Iopt [struct]: Structure of identification options. See idModels.NsfPolyModel.identifyModelOptions()

% Data
T_out = Data.WeatherData.data(:,2); % Outside Temp 
P_sun = Data.WeatherData.data(:,9)/1e3; % Global hor rad. [kWh/m2]

% P.GeoPos = [52.47 13.40]; 
% H = Hvaccalc.calcSunPosition((0:(length(T_out)-1))/24+datenum('2020-1-1'),P.GeoPos(2),P.GeoPos(1));
% H = (H>0).*H;


% Model Setup
Mdl.Ts = 3600;
Mdl.TimeUnit = 'seconds';
if ~isempty(P.fc)
    Mdl.PreFilter = cell(1,2); 
    [Mdl.PreFilter{:}] = butter(1,(P.fc/60)/((1/(Mdl.Ts/60))/2),'high');
    Mdl.PreFilter{1} = Mdl.PreFilter{1}./Mdl.PreFilter{1}(1);
end
Mdl.OutputName = {'Außentemperatur'; 'GlobalStrahlung'};
Mdl = util.struct2namevaluepairs(Mdl);

% Do identification (Outside temperature)
% Mdl.OutputName  = 'Außentemperatur';
% [R,M{1}] = idModels.test.sweepPar({'Order'},@idModels.NsfPolyModel.estModel,T_out,u,'Order',P.Orders,...
%                             'Type',P.Type,'InitFromPreviousN',P.InitFromPrevious,'MdlProp',Mdl,...
%                             'HpVal',P.HpVal,'FactorizeG',P.FactorizeG,'IdOpt',P.Iopt,...
%                             'SeasonalityC',P.Seasonality,'SeasonalityD',P.Seasonality,...
%                             'SeasonalityOrder',P.SeasonalityOrder,...
%                             'MinRealTol',P.MinRealTol,'ObserverInitSamples','auto');
% [mat(:,1), lab] = tab2mat(R,1);    
% [~, ix_srt] = sort(R.Order);
% M{1} = M{1}(ix_srt);
% 
% % Do identification Radiation
% Mdl.OutputName  = 'Globalstrahlung';
% [R,M{2}] = idModels.test.sweepPar({'Order'},@idModels.NsfPolyModel.estModel,P_sun,u,'Order',P.Orders,...
%                             'Type',P.Type,'InitFromPreviousN',P.InitFromPrevious,'MdlProp',Mdl,...
%                             'HpVal',P.HpVal,'FactorizeG',P.FactorizeG,'IdOpt',P.Iopt,...
%                             'SeasonalityC',P.Seasonality,'SeasonalityD',P.Seasonality,...
%                             'SeasonalityOrder',P.SeasonalityOrder,...
%                             'MinRealTol',P.MinRealTol,'ObserverInitSamples','auto');
% [mat(:,2), lab] = tab2mat(R,1);    
% [~, ix_srt] = sort(R.Order);
% M{2} = M{2}(ix_srt);
% 
% % Order selection and output models
% figure('color','w','Position',[200 200 600 300]); 
% subplot(1,2,1); bar(mat(:,1)); xticklabels(arrayfun(@num2str,P.Orders,'UniformOutput',0)); 
% subplot(1,2,2); bar(mat(:,2)); xticklabels(arrayfun(@num2str,P.Orders,'UniformOutput',0)); 
% ylab = {'RMMSE ($\varepsilon [t|t-k]$) [K]' 'RMMSE ($\varepsilon [t|t-k]$) [kW/m$^2$]'};
% formatFigure(14,'Polynomordnung',ylab,[],{'Au{\ss}entemp.' 'Globalstrahlung'});
% 
% ix = input('Choose Model order for outside temperature model!) \n');
% Mdl_Tout = M{1}{ix};
% ix = input('Choose Model order for radiation model!) \n');
% Mdl_Psun = M{2}{ix};

% ix = input('Choose Model order for outside temperature model!) \n');
% Mdl_Tout = M{1}{ix};
% ix = input('Choose Model order for radiation model!) \n');
% Mdl_Psun = M{2}{ix};

if ischar(P.SeasonalityOrder) && strcmpi(P.SeasonalityOrder,'ByOrder')
	Inc = true;
else
	Inc = false;
end

for k = 1:length(P.Orders)
    if ~isempty(P.Seasonality)
        ix = [false true(1,P.Seasonality+Inc*(P.Orders(k)-1))];
        ix(P.Orders(k)+2:P.Seasonality-Inc*(P.Orders(k)-1)) = false;
    else
        ix = [false true(1,P.Orders(k))];
    end
    nA = (length(ix)-1)*ones(2);
    nA(2,1) = 0;
    nA(1,2) = 0;
    M{k} = idModels.NsfPolyModel(nA,[],[],(length(ix)-1)*ones(2,1),Mdl{:});
    if ~isempty(P.Seasonality)
        for nr = 1:size(M{k}.A,1)
            M{k}.C(nr).val([false ~ix(2:end)]) = 0;
            M{k}.C(nr).free = ix;
            %for nc = 1:size(M{k}.A,2)
              	M{k}.A(nr,nr).val([false ~ix(2:end)]) = 0;
                M{k}.A(nr,nr).free = ix;
            %end
        end
    end
    P.Iopt.Hp = 1;
    M{k}.identify([T_out P_sun],[],P.Iopt);
    M{k}.updateSs('minrealtol',1e-8);
    e = M{k}.calcResiduals([T_out P_sun],[],'Hp',P.HpVal);
    for no = 1:size(e,2)
        E = squeeze(e(:,no,:));
        mat(k,no) = mean(abs(E(:)),'omitnan');
    end
end

% Order selection and output models
figure('color','w','Position',[200 200 600 300]); 
subplot(1,2,1); bar(mat(:,1)); xticklabels(arrayfun(@num2str,P.Orders,'UniformOutput',0)); 
subplot(1,2,2); bar(mat(:,2)); xticklabels(arrayfun(@num2str,P.Orders,'UniformOutput',0)); 
ylab = {'RMMSE ($\varepsilon [t|t-k]$) [K]' 'RMMSE ($\varepsilon [t|t-k]$) [kW/m$^2$]'};
util.formatFigure(14,'Polynomordnung',ylab,[],{'Au{\ss}entemp.' 'Globalstrahlung'});
idx = input('Choose Model order!) \n');
M = M{idx};
end
