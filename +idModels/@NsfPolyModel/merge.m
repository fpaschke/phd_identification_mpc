function M = merge(m1,varargin)
%MERGE Merges multiple Models to one model. If an input of Model i is an 
%output of model j the function tries to merge the Model by replacing the inputs by
%an corresponding outputfeedback. Since the implementation of the Inputnonlienarity 
%function handles is not known the user will be prompted for an input.
%(exception: if 'NewOutputIdx' is supplied then user will not be prompted).
%
% M = merge(m1,m2,...,mN,varargin)
% M [idModels.NsfPolyModel]: Merged model
% m1,m2,...,mN [idModels.NsfPolyModel]: The models to be merged
% varargin: optional Name Value Pairs:
%   'Name' [char]: The Name of the new model. Def. 'unnamed' 
%   'NewOutputIdx' [cell of int vectors]:   Only necessary when input of Model i is an output of model j. Def. {}
%                                           The user will not be prompted if NewOutputIdx is supplied. Instead 
%                                           The Parametervalues will be used as new output_idx.

%% Parse Inputs
ix_firstchar = find(cellfun(@ischar,varargin),1,'first');
p = inputParser();
addParameter(p,'Name','unnamed',@ischar)
addParameter(p,'NewOutputIdx',{},@iscell)
parse(p,varargin{ix_firstchar:end});
if ~isempty(ix_firstchar); m2add = ix_firstchar-1; else; m2add = length(varargin); end % Get number of Models to be merged Models

%% Init new Properties
tu = m1.TimeUnit;
In_new = m1.InputName;
Out_new = m1.OutputName;
G_new = m1.InputNonlinearity; % the new Manipulated InputNonlinearity structure
A_new = m1.A;
B_new = m1.B;
C_new = m1.C;
D_new = m1.D;
E_new = m1.E;
F_new = m1.F;
Pf_new = m1.PreFilter;
y0_new = m1.OutputOffset;
Nv_new = m1.NoiseVariance;
InMax_new = m1.InputMax;
InMin_new = m1.InputMin;
OutMin_new = m1.OutputMin;
OutMax_new = m1.OutputMax;

%% Start Merging
for mi = 1:m2add
    Mi = varargin{mi};
	assert(Mi.OutputDimension == 1,'Merging of MIMO Models not supported yet!');
 	assert(isequal(Mi.TimeUnit,tu),'TimeUnits of models are not equal!');
	assert(m1.Ts==Mi.Ts,'Models do not have equal sampling times. Consider resampling first!');
	assert(isequal(length(unique(Out_new)),length(Out_new)),'The Names of the Outputs need to be different!');
    ny = length(Out_new)+1:length(Out_new)+Mi.OutputDimension; % the indexes of the new outputs wrt the new model
	Out_new = [Out_new; Mi.OutputName];
    y0_new = [y0_new; Mi.OutputOffset];
    A_new(ny,ny) = Mi.A; % nebendiagonalen!!!
    C_new = [C_new; Mi.C];
    D_new = [D_new; Mi.D];
    E_new = [E_new; Mi.E];
    Pf_new = [Pf_new; Mi.PreFilter];
    Nv_new(ny,ny) = Mi.NoiseVariance;
    OutMin_new = [OutMin_new; Mi.OutputMin];
    OutMax_new = [OutMax_new; Mi.OutputMax];
    
    % Now append Mi.InputNonlinearity to G_new, Mi.B and Mi.F
    for k = 1:length(Mi.InputNonlinearity)
        g_k = Mi.InputNonlinearity(k);
        if ~isempty(g_k.output_idx) %we need to exchange the output_idxs
            for no = 1:length(g_k.output_idx)
                g_k.output_idx(no) = find(strcmp(Mi.OutputName{g_k.output_idx(no)},Out_new)); % the inputnumber of Mi.InputName{ni} wrt In_new 
            end
        end
        for ni = 1:length(g_k.input_idx)
            thisname = Mi.InputName{g_k.input_idx(ni)};
            ix_in = find(strcmp(thisname,In_new)); % the inputnumber of Mi.InputName{ni} wrt In_new 
            if ~isempty(ix_in) 
                InMin_new(ix_in) = min(InMin_new(ix_in),Mi.InputMin(g_k.input_idx(ni)));
                InMax_new(ix_in) = max(InMax_new(ix_in),Mi.InputMax(g_k.input_idx(ni)));
                g_k.input_idx(ni) = ix_in; %this input is defined already -> change index
            else % Mi obviously has additional input -> add
                In_new = [In_new; thisname]; % add new inputnames and maxima/minima
                InMax_new = [InMax_new; Mi.InputMax(g_k.input_idx(ni))];
                InMin_new = [InMin_new; Mi.InputMin(g_k.input_idx(ni))];
                g_k.input_idx(ni) = length(In_new); % replace by the new index
            end
        end

        % Now check weather g_k matches somebody in G_new - > we do not need to add in this case
        ix_equal = 0;
        for j = 1:length(G_new)
            ix_equal = isequal(G_new(j),g_k);
            if ix_equal % ok we need to put B and F to the correct place
                B_new(ny,j) = Mi.B(:,k); 
                F_new(ny,j) = Mi.F(:,k); 
                break;
            end
        end
        if ~ix_equal % ok g_k isnt contained
            G_new = [G_new; g_k];
            B_new(ny,size(B_new,2)+1) = Mi.B(:,k); 
            F_new(ny,size(F_new,2)+1) = Mi.F(:,k); 
        end
    end
end

%% Print Output channels
fprintf('NEW MODEL OUTPUTS: \n')
for no = 1:length(Out_new)
    fprintf('Output Channel %i:\t %s \n',no,Out_new{no}); 
end

%% Now check weather inputnames Contains Outputnames
k = 1;
for no = 1:length(Out_new)
    ix = find(strcmp(In_new,Out_new{no})); % InputName{ix} is actually OutputName{no}!
    if ~isempty(ix) % ok we need to too iterate through all inputnonlinearities and look for input_idx ix
      	% print current inputs
        if isempty(p.Results.NewOutputIdx) 
            fprintf('CURRENT MODEL INPUTS: \n')
            for ni = 1:length(In_new)
                fprintf('Input Channel %i:\t %s \n',ni,In_new{ni}); 
            end
        end
        % iteratre through all iputnonlineartities
        for ni = 1:length(G_new) 
            if any(G_new(ni).input_idx==ix)
                % We need to ask the user how output_idx should look like after merging the models
                old_length = length(G_new(ni).output_idx);
                if isempty(p.Results.NewOutputIdx) 
                    fprintf('User input is necessary. The Inputnonlinearity %i: \n',ni);
                    G_new(ni)
                    fprintf('contains inputchannel %i (%s) which actually is output number %i (%s) of the merged model. \n',ix,In_new{ix},no,Out_new{no});
                    fprintf('This input will be deleted from input_idx and output_idx will be modified! \n');
                    fprintf('However it is not clear how output_idx should be arraged since it depends on the implementation of the InputNonlinearity function handle! \n');
                    fprintf('Please enter the new vector output_idx (It should have length %i):\n',old_length+1);
                end
                while 1
                    if ~isempty(p.Results.NewOutputIdx) 
                        G_new(ni).output_idx = p.Results.NewOutputIdx{k};
                    else
                        G_new(ni).output_idx = input('output_idx = ');
                    end
                    if all(G_new(ni).output_idx<=length(Out_new)) && length(G_new(ni).output_idx) == old_length + 1
                        break;
                    else
                        frintf('ERROR: Assignment is invalid! The new output_idx should have length %i and only values between 1 and %i',old_length + 1,length(Out_new));
                    end
                end
                % Delete ix from input_idx
                G_new(ni).input_idx(G_new(ni).input_idx==ix) = [];
                % increase counter
                k = k + 1;
            end
        	% All input indexes > ix need to be decreased by 1
            G_new(ni).input_idx(G_new(ni).input_idx>ix) = G_new(ni).input_idx(G_new(ni).input_idx>ix) - 1;
        end
        % Now remove InputName
        In_new(ix) = [];
        InMax_new(ix) = [];
        InMin_new(ix) = [];
    end
end

%% Print new Model inputs
fprintf('NEW MODEL INPUTS: \n')
for ni = 1:length(In_new)
    fprintf('Input Channel %i:\t %s \n',ni,In_new{ni}); 
end

%% Now Initialize new Model
M = idModels.NsfPolyModel(zeros(size(A_new)),zeros(size(B_new)),zeros(size(B_new)),zeros(size(C_new)),zeros(size(D_new)),zeros(size(E_new)),zeros(size(F_new)),...
                                'InputNonlinearity',G_new,'Ts',m1.Ts,'TimeUnit',m1.TimeUnit,'InputName',In_new,'OutputName',Out_new,...
                                'A',A_new,'B',B_new,'C',C_new,'D',D_new,'E',E_new,'F',F_new,...
                                'OutputOffset',y0_new,'Name',p.Results.Name,'PreFilter',Pf_new,...
                                'NoiseVariance',Nv_new,'InputMin',InMin_new,'InputMax',InMax_new,'OutputMax',OutMax_new,'OutputMin',OutMin_new);
end
