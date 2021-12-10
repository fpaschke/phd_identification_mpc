function Mdl_out = estRoomModel(Data,P)
% ESTROOMMODEL Setups the Model structure, estimates models from data
% and displays the identification result. At the end user will be prompted
% to select model orders. See also: idModels.NsfPolyModel.estModel for more
% information on parameters.
%
% Mdl_out [4 x 1 cell of idModels.NsfPolyModel]: Estimated models
% Data [struct or object]:  Contains timeseries data necessary for identification (N #Samples):
%                           All Signals are assumed to be uniformly sampled.
%           Data.T_air [4 x N]: Air Temperatures for each of the rooms [C] (output)
%           Data.T_sup [1 x N]: Supply temp of AHU [C] (input)
%           Data.m_sup [4 x N]: Supply air massflowrates to each room [kg/s] (input)
%           Data.T_out [1 x N]: Outside Air Temperature [C] (input)
%           Data.P_sun [1 x N]: Global Radiation on horizontal plane [kW/m2] (input)
%           Data.Azimuth [1 x N]: Azimuth of sun [degrees] (input)
%           Data.Elevation [1 x N]: Height of sun [degrees] (input)
% P [struct]: Structure containing Model properties validation and identification options.
%           P.WinDir [4 x 1 double]:Orientation of windows wrt to North (here: [270 90 270 90])
%           P.fc [double]:          Cutoff Frequency of 1st order highpass butterworth prefilter in [1/h]. 
%                                   1/(6*2*pi) corresponds to 1st Order System with Timeconstant of 6h
%           P.Orders [int vector]:  Orders of Polynomials A,B,... for which identification will be performed
%           P.Type [char]:          Model Structure that will be used for identification ('arx', 'oe', 'armax' or 'bj')  
%           P.FactorizeG [logical]:	Identifiy model with factorized A(q),B(q),E(q) if true
%           P.MinRealTol [logical]: Perform model order reduction of ss object using minreal(ss,tol)
%           P.Seasonality (int or []): Seasonality of NoiseModel in number of samples
%           P.SeasonalityOrder [logical]: Seasonality will be Increased by Order if 'ByOrder'
%           P.HpVal [int]:          Prediction horizon in #steps for validation (plot)
%           P.InitFromPrevious [log]: Use previous estimated Model for Parameter initialization
%           P.Iopt [struct]: Structure of identification options. See idModels.NsfPolyModel.identifyModelOptions()

% Model Setup
iopt = util.struct2namevaluepairs(P.Iopt);
Nz = size(Data.m_sup.Data,2);
Mdl.Ts =  P.T_s*60;
Mdl.TimeUnit = Data.T_sup.TimeInfo.Units;
assert(strcmp(Mdl.TimeUnit,'seconds'),'Timeunits of timeseries need to be in seconds!');

% Input Vector
N = length(Data.T_sup.Data);
dk = Mdl.Ts/Data.T_air.TimeInfo.Increment;
u_ = [Data.In.T_out.Data(1:dk:N) ...
      Data.In.P_sun.Data(1:dk:N) ...
      Data.In.Azimuth.Data(1:dk:N) ...
      Data.In.Elevation.Data(1:dk:N) ...
      Data.T_sup.Data(1:dk:N)];

ny = size(P.Orders{1},1);
Ord = [P.Orders mat2cell(zeros(ny,6-length(P.Orders)),ny,ones(6-length(P.Orders),1))];
if ~isempty(P.Seasonality)
    for no = 1:max(1,ny-1)
        ixC{no} = [false true(1,P.Seasonality+P.SeasonalityOrder-1)];
        ixC{no}(Ord{4}(no)+2:P.Seasonality-P.SeasonalityOrder+1) = false;
        Ord{4}(no) = (P.Seasonality+P.SeasonalityOrder-1);
        ixD{no} = [false true(1,P.Seasonality+P.SeasonalityOrder-1)];
        ixD{no}(Ord{5}(no)+2:P.Seasonality-P.SeasonalityOrder+1) = false;
        Ord{5}(no) = (P.Seasonality+P.SeasonalityOrder-1);
    end
end
  
if P.MIMO
    u = [u_ Data.m_sup.Data(1:dk:end,:)];
    y = Data.T_air.Data(1:dk:end,:);
    Mdl.InputName = {'Außentemperatur' 'Globalstrahlung' 'Azimuth' 'Sonnenhöhe' 'Zulufttemperatur' 'Massenstrom HS1' 'Massenstrom HS2' 'Massenstrom HS3' 'Massenstrom HS4'}; 
    Mdl.OutputName  = {'Lufttemperatur HS1' 'Lufttemperatur HS2' 'Lufttemperatur HS3' 'Lufttemperatur HS4' 'Lufttemperatur Flur'};
    for nz = 1:Nz+1
        Mdl.InputNonlinearity(nz+Nz).fun = 'idModels.func.fun_rad'; 	% calculates radiation on inclined plane 
        Mdl.InputNonlinearity(nz+Nz).input_idx = [2 3 4];               % Radiation, Azimuth, Height 
        Mdl.InputNonlinearity(nz+Nz).free = false(1,3); 
        if nz <= Nz
            Mdl.InputNonlinearity(nz+Nz).parameters = [P.WinDir(nz) 90 2];
            Mdl.InputNonlinearity(nz).fun = 'idModels.func.fun_VdT';	% virt heating Power
            Mdl.InputNonlinearity(nz).parameters = [1 0];               % m*(Tsup-p(2)-Troom)^p(1)
            Mdl.InputNonlinearity(nz).free = false(1,2);              
            Mdl.InputNonlinearity(nz).input_idx = [5+nz 5];          	% Vol Supply Temp   
            Mdl.InputNonlinearity(nz).output_idx = nz;                  % Room temp    
        else
            Mdl.InputNonlinearity(nz+Nz) = Mdl.InputNonlinearity(nz+Nz-1);      % Northern entrence
         	Mdl.InputNonlinearity(nz+Nz+1) = Mdl.InputNonlinearity(nz+Nz-1);	% Southern entrence
            Mdl.InputNonlinearity(nz+Nz).parameters = [0 90 2];         
            Mdl.InputNonlinearity(nz+Nz+1).parameters = [180 90 2];        
        end
    end
    Mdl.InputNonlinearity(11).fun = [];                                 % Identity
    Mdl.InputNonlinearity(11).input_idx = 1;                            % Outside Temp     
    
    mdl_opt = util.struct2namevaluepairs(Mdl);
    if iscell(P.Orders)
        Mdl_out = idModels.NsfPolyModel(Ord{:},mdl_opt{:});
        if ~isempty(P.Seasonality)
            for no = 1:ny-1
                Mdl_out.C(no).val([false ~ixC{no}(2:end)]) = 0;
                Mdl_out.D(no).val([false ~ixD{no}(2:end)]) = 0;
                Mdl_out.C(no).free = ixC{no};
                Mdl_out.D(no).free = ixD{no};
            end
        end
        if P.FactorizeG
            Mdl_out.factorize({'A' 'B' 'E'});
        end
        Mdl_out.identify(y,u,iopt{:});
     	Mdl_out.updateSs('BalRedOrder',30); %% P/Z Cancellation using (minreal)
    	if ~isempty(P.Seasonality)
            Mdl_out.updateSs('Sparse',true)
        end
        e = Mdl_out.calcResiduals(y,u,'Hp',P.HpVal); % InitMethod of observer = Inf; no GLS if true 
        for nz = 1:Nz+1
            e_nz = squeeze(e(:,nz,:));
            e_nz = e_nz(:);
            E(1,nz) = mean(abs(e_nz),'omitnan');
            fprintf('RMMSE Hs %i: %f \n',nz,E(nz));
        end
    else
        for k = 1:length(P.Orders)        
            n = P.Orders(k);
            nA =   n*[  1 0 1 0 1; 
                        0 1 0 1 1; 
                        1 0 1 0 1;
                        0 1 0 1 1;
                        1 1 1 1 1];

            nB =   n*[  1 0 0 0 1 0 0 0 0 0 1; 
                        0 1 0 0 0 1 0 0 0 0 1; 
                        0 0 1 0 0 0 1 0 0 0 1;
                        0 0 0 1 0 0 0 1 0 0 1;
                        0 0 0 0 0 0 0 0 1 1 1];
            M{k} = idModels.NsfPolyModel(nA,nB,ones(size(nB)),n*ones(Nz+1,1),mdl_opt{:});
            M{k}.identify(y,u,iopt{:});
            if ~isempty(P.Seasonality)
                M{k}.updateSs('Sparse',true)
            end
            e = M{k}.calcResiduals(y,u,'Hp',P.HpVal); % InitMethod of observer = Inf; no GLS if true 
            for nz = 1:Nz+1
                e_nz = squeeze(e(:,nz,:));
                e_nz = e_nz(:);
                E(n,nz) = mean(abs(e_nz),'omitnan');
            end
        end
        figure('color','w'); bar(E); xticklabels(arrayfun(@num2str,P.Orders,'UniformOutput',0));
        formatFigure(14,'Polynomordnung','RMMSE ($\varepsilon [t|t-k]$)'); legend({'HS 1' 'HS2' 'HS3' 'HS4' 'Flur'});
        ix = input('Choose Model order \n');
        Mdl_out = M{ix};
    end
else
    Mdl.InputNonlinearity(1).fun = [];                          % Identity
    Mdl.InputNonlinearity(1).input_idx = 1;                 	% Outside Temp        
    Mdl.InputNonlinearity(2).fun = 'idModels.func.fun_rad';     % calculates radiation on inclined plane 
    Mdl.InputNonlinearity(2).input_idx = [2 3 4];               % Radiation, Azimuth, Height
    if any(any(Data.m_sup.Data>0))
        Mdl.InputNonlinearity(3).fun = 'idModels.func.fun_VdT';	% virt heating Power
        Mdl.InputNonlinearity(3).parameters = [1 0];          	% m*(Tsup-p(2)-Troom)^p(1)
        Mdl.InputNonlinearity(3).free = [false false];              
        Mdl.InputNonlinearity(3).input_idx = [6 5];          	% Vol Supply Temp   
        Mdl.InputNonlinearity(3).output_idx = [1];              
    end
    % Prefilter (1st order butterworth)
    if ~isempty(P.fc)
        Mdl.PreFilter = cell(1,2); 
        [Mdl.PreFilter{:}] = butter(1,(P.fc/60)/((1/(Mdl.Ts/60))/2),'high');
        Mdl.PreFilter{1} = Mdl.PreFilter{1}./Mdl.PreFilter{1}(1);
    end
    for nhs = 1:Nz
        % Names
        Mdl.InputName = {'Außentemperatur' 'Globalstrahlung' 'Azimuth' 'Sonnenhöhe' 'Zulufttemperatur' ['Massenstrom HS' num2str(nhs)]}; 
        Mdl.OutputName  = ['Lufttemperatur HS' num2str(nhs)];

        % orientation of windows wrt N in degrees, orientation of windows wrt to horizontal plane in degrees
        Mdl.InputNonlinearity(2).parameters = [P.WinDir(nhs) 90 2];     
        Mdl.InputNonlinearity(2).free = false(1,3); 

        % Build input and output vector
        y = Data.T_air.Data(1:dk:N,nhs);
        u = [u_ Data.m_sup.Data(1:dk:N,nhs)];
        if ~iscell(P.Orders) % Try different orders
            % Do identification
            [R{nhs},M{nhs}] = idModels.test.sweepPar({'Order'},@idModels.NsfPolyModel.estModel,y,u,'Order',P.Orders,...
                                        'Type',P.Type,'InitFromPreviousN',P.InitFromPrevious,'MdlProp',Mdl,...
                                        'HpVal',P.HpVal,'FactorizeG',P.FactorizeG,'IdOpt',P.Iopt,...
                                        'SeasonalityC',P.Seasonality,'SeasonalityD',P.Seasonality,...
                                        'SeasonalityOrder',P.SeasonalityOrder,...
                                        'MinRealTol',P.MinRealTol,'ObserverInitSamples','auto');
            [mat(:,nhs), lab] = tab2mat(R{nhs},1);    
            [~, ix_srt] = sort(R{nhs}.Order);
            M{nhs} = M{nhs}(ix_srt);
            leg{nhs} = ['Room ' num2str(nhs)];
        else
            mopt = util.struct2namevaluepairs(Mdl); 
            M{nhs} = idModels.NsfPolyModel(Ord{:},mopt{:});
            if ~isempty(P.Seasonality)
                M{nhs}.C.val([false ~ixC{1}(2:end)]) = 0;
                M{nhs}.D.val([false ~ixD{1}(2:end)]) = 0;
                M{nhs}.C.free = ixC{1};
                M{nhs}.D.free = ixD{1};
            end
            if P.FactorizeG
                M{nhs}.factorize({'A' 'B' 'E'});
            end
            M{nhs}.identify(y,u,iopt{:});
            e = M{nhs}.calcResiduals(y,u,P.HpVal); % InitMethod of observer = Inf; no GLS if true 
            RMMSE(nhs) = sqrt(mean(e(:).^2,'omitnan'));
        end
    end

    % Order selection and output models
    if ~iscell(P.Orders)
        figure('color','w'); bar(mat); xticklabels(arrayfun(@num2str,P.Orders,'UniformOutput',0));
        formatFigure(14,'Polynomordnung','RMMSE ($\varepsilon [t|t-k]$)',[],[],leg);
        ix = input('Choose Model order (Enter 4 element vector of intergers corresponding to model order to be chosen for each hall!) \n');
    end
    for nhs = 1:Nz
        if ~iscell(P.Orders)
            M{nhs} = M{nhs}{ix(nhs)}; 
        else
            fprintf('RMMSE Hs %i: %f K\n',nhs,RMMSE(nhs));
        end
        if P.FactorizeG
            if contains(M{nhs}.Type{1},'AR') 
                rts = M{nhs}.A.val(2:end);
            else
                rts = M{nhs}.E.val(2:end);
            end
            fprintf('TimeConstant(s) Room %i: %s h\n',nhs,mat2str(-Mdl.Ts./log(rts)/3600,4));
        end
    end
    Mdl_out = merge(M{:});
    Mdl_out.Info.OptionsUsed.MultiStep = M{nhs}.Info.OptionsUsed.MultiStep;
    Mdl_out.Info.OptionsUsed.Hp = M{nhs}.Info.OptionsUsed.Hp;
end
