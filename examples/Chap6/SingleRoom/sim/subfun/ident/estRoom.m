function Mdl_out = estRoom(Data,P)
% ESTROOM Setups the Model structure, estimates models from data
% and displays the identification result. At the end user will be prompted
% to select model order. See also: idModels.NsfPolyModel.estModel for more
% information on parameters.
%
% Mdl_out [4 x 1 cell of idModels.NsfPolyModel]: Estimated models
% Data [SimulationOutput]:  Contains timeseries data necessary for identification.
%                           Will be generated by simulate()
% P [struct]: Structure containing Model properties validation and identification options.
%           P.WinDir [double]:  Orientation of window wrt to North (here: [270])
%           P.fc [double]:          Cutoff Frequency of 1st order highpass butterworth prefilter in [1/h]. 
%                                   1/(6*2*pi) corresponds to 1st Order System with Timeconstant of 6h
%           P.Orders [int vector]:  Orders of Polynomials A,B,... for which identification will be performed
%           P.Type [char]:          Model Structure that will be used for identification ('arx', 'oe', 'armax' or 'bj')  
%           P.FactorizeG [logical]:	Identifiy model with factorized A(q),B(q),E(q) if true
%           P.HpVal [int]:          Prediction horizon in #steps for validation (plot)
%           P.InitFromPrevious [log]: Use previous estimated Model for Parameter initialization
%           P.Iopt [struct]: Structure of identification options. See idModels.NsfPolyModel.identifyModelOptions()

%% Identify Room Model 
Mdl = struct;
Mdl.Name = 'Room Model';
Mdl.Ts =  Data.T_air.TimeInfo.Increment;
Mdl.TimeUnit = 'seconds';

Mdl.InputName = {   Data.In.T_out.Name 
                    Data.In.P_sun.Name 
                    Data.In.Azimuth.Name 
                    Data.In.Elevation.Name 
                    Data.H_win.Name
                    Data.H_flo.Name
                    'Vorlauftemperatur'
                    %Data.In.P_sun_diff.Name
                    }; 
Mdl.OutputName  = 'Raumtemperatur';
Mdl.InputNonlinearity(1).fun = ''; % Identity function for outside temp.
Mdl.InputNonlinearity(1).parameters = [];
Mdl.InputNonlinearity(1).free = []; 
Mdl.InputNonlinearity(1).input_idx = 1; % OutsideTemp
Mdl.InputNonlinearity(1).output_idx = [];
Mdl.InputNonlinearity(2).fun = 'idModels.func.fun_rad'; % calculates radiation on inclined plane 
Mdl.InputNonlinearity(2).parameters = [P.WinDir 90 2]; % orientation of windows wrt N in degrees, orientation of windows wrt to horizontal plane in degrees, method to compute diffuse radiation 
Mdl.InputNonlinearity(2).free = [false false false];
Mdl.InputNonlinearity(2).input_idx = [2 3 4 5]; % Radiation, Azimuth, Height, Blind pos   
Mdl.InputNonlinearity(2).output_idx = [];
Mdl.InputNonlinearity(3).fun = 'idModels.func.fun_powdT'; % "virt. Heatingpower": mean(Valvepos)^p1  * (ply_heat - Troom)^p2 
Mdl.InputNonlinearity(3).parameters = [1 1.1]; % initial estimates of p1 and p2
Mdl.InputNonlinearity(3).free = [true false]; 
Mdl.InputNonlinearity(3).input_idx = [6 7]; % Valvepos, SupplyTemp
Mdl.InputNonlinearity(3).output_idx = 1;
if false % Only for testing
    Mdl.InputNonlinearity(4).fun = 'idModels.func.fun_rad_diff'; % calculates radiation on inclined plane and multiply by logistic function to consider shading by trees and buildings in front of the window
    Mdl.InputNonlinearity(4).parameters = []; % orientation of windows wrt N in degrees, orientation of windows wrt to horizontal plane in degrees, estimated shading angle, growth rate/steepness of logistic function
    Mdl.InputNonlinearity(4).free = []; 
    Mdl.InputNonlinearity(4).input_idx = [2 5]; % Radiation, Blind pos   
    Mdl.InputNonlinearity(4).output_idx = [];
end

if ~isempty(P.fc)
    Mdl.PreFilter = cell(1,2); 
    if (P.fc/60)/((1/(Mdl.Ts/60))/2)>1
        Mdl.PreFilter = {[1 -1] 1};
    else
        [Mdl.PreFilter{:}] = butter(1,(P.fc/60)/((1/(Mdl.Ts/60))/2),'high');
        Mdl.PreFilter{1} = Mdl.PreFilter{1}./Mdl.PreFilter{1}(1);
    end
end
ns = 1;
N = length(Data.T_sup_flo.Data);
u = [Data.In.T_out.Data(ns:N) ...
    Data.In.P_sun.Data(ns:N) ...
    Data.In.Azimuth.Data(ns:N) ...
    Data.In.Elevation.Data(ns:N)...
    Data.H_win.Data(ns:N) ... 
    Data.H_flo.Data(ns:N) ... 
    Data.T_sup_flo.Data(ns:end) ...
    %Data.In.P_sun_diff.Data(ns:N) ...
    ];
y = Data.T_air.Data(ns:end);

CostType = 1;
if ~isfield(P,'Seasonality'); P.Seasonality = []; end
if ~isfield(P,'SeasonalityOrder'); P.SeasonalityOrder = 0; end
if ~iscell(P.Orders) % Try different orders
    varpar = {'Order'};
    if iscell(P.Type) && length(P.Type)>1
        varpar = [varpar 'Type'];
    end
    [R,M] = idModels.test.sweepPar(varpar,@idModels.NsfPolyModel.estModel,y,u,'Order',P.Orders,...
                'Type',P.Type,'InitFromPreviousN',P.InitFromPrevious,'MdlProp',Mdl,...
                'HpVal',P.HpVal,'FactorizeG',P.FactorizeG,'IdOpt',P.Iopt,'ObserverInitSamples','auto',...
                'SeasonalityC',P.Seasonality,'SeasonalityD',P.Seasonality);
    [mat, lab] = tab2mat(R,CostType);    
    [~, ix_srt] = sort(R.Order);
    M = M(ix_srt);

    % Order selection and output models
    figure('color','w','Position',[200 200 350 280]); bar(mat); xticklabels(arrayfun(@num2str,P.Orders,'UniformOutput',0));
    formatFigure(14,'Polynomordnung','RMMSE ($\varepsilon [t|t-k]$)');
    if iscell(P.Type) && length(P.Type)>1
        legend(upper(P.Type),'Interpreter','latex');
    end
    ix = input('Choose Model order! \n');
    Mdl_out = M{ix}; 
else % Use Specified order
    mopt = util.struct2namevaluepairs(Mdl); 
    iopt = util.struct2namevaluepairs(P.Iopt);
	if ~isempty(P.Seasonality)
        nC = P.Orders{4}; nD = P.Orders{5};
        if strcmpi(P.SeasonalityOrder,'ByOrder') 
            if P.Orders{1}~=0
                nS = P.Orders{1}-1;
            else
                nS = P.Orders{6}-1;
            end
        else
            nS = 0;
        end
        P.Orders{4} = P.Seasonality+nS; %C
        P.Orders{5} = P.Seasonality+nS; %D
  	end
    Mdl_out = idModels.NsfPolyModel(P.Orders{:},mopt{:});
    if ~isempty(P.Seasonality)
      	Mdl_out.C.val(nC+2:P.Seasonality-nS) = 0;
        Mdl_out.C.free(nC+2:P.Seasonality-nS) = false;
     	Mdl_out.D.val(nD+2:P.Seasonality-nS) = 0;
        Mdl_out.D.free(nD+2:P.Seasonality-nS) = false;
    end
  	if P.FactorizeG
        Mdl_out.factorize({'A' 'B' 'E'});
    end
    Mdl_out.identify(y,u,iopt{:});
 	e = Mdl_out.calcResiduals(y,u,P.HpVal); % InitMethod of observer = Inf; no GLS if true 
    if CostType == 1
        fprintf('RMMSE: %f K\n',sqrt(mean(e(:).^2,'omitnan')));
    else
        fprintf('RMMSE: %f K\n',sqrt(mean(e(:,:,end).^2,'omitnan')));
    end
end

if P.FactorizeG
    if contains(Mdl_out.Type,'ARMAX') || contains(Mdl_out.Type,'ARX') || contains(Mdl_out.Type,'ARIX') || contains(Mdl_out.Type,'ARIMAX')
        rts = Mdl_out.A.val(2:end);
    else
        rts = Mdl_out.E.val(2:end);
    end
end
fprintf('TimeConstant Room: %s min\n',mat2str(-Mdl.Ts./log(rts)/60,4));
end
