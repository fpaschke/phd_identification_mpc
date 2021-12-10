function MDL = merge(varargin)
%MERGE Merges multiple GeneralizedStaticMopdels into one model. 
%Only Models with equal inputs i.e. obj.InputName property can be merged.
%
% MDL = merge(mdl1,mdl2,...,mdlN)
% mdl1 ... mdlN [idModels.GeneralizedStaticModel]: Models to be merged.

Nm = length(varargin);
InNames = varargin{1}.InputName; OutNames = {}; Params = {}; Free = {}; Nv = []; InMax = []; InMin = []; OutMax = []; OutMin = []; 
OutOff = []; Fun = {}; Cl = {};
for nm = 1:Nm
	assert(isequal(varargin{nm}.InputName,InNames),'Models can be Merged only if they have equal InputName!');
    assert(isa(varargin{nm},'idModels.GeneralizedStaticModel'),'All models need to be GeneralizedStaticModels!');
    assert(varargin{nm}.OutputDimension == 1,'Currently only one output models can be merged!');
    Cl = [Cl; class(varargin{nm})];
    Fun = [Fun; varargin{nm}.Fun];
    OutNames = [OutNames; varargin{nm}.OutputName];
    Nv = [Nv; varargin{nm}.NoiseVariance];
    Params = [Params; varargin{nm}.Parameters];
    Free = [Free; varargin{nm}.Free];
    InMax = [InMax varargin{nm}.InputMax];
    InMin = [InMin varargin{nm}.InputMin];
    OutMax = [OutMax; varargin{nm}.OutputMax];
    OutMin = [OutMin; varargin{nm}.OutputMin];
    OutOff = [OutOff; varargin{nm}.OutputOffset];
end
assert(length(varargin) == length(unique(OutNames)),'All OutputNames need to be different!');
if length(unique(Cl))==1
	MDL = eval([Cl{1} '(''InputName'',varargin{1}.InputName,''OutputName'',OutNames)']);
    % HANDLING OF SUBCLASS PROPERIES!
else
    MDL = idModels.GeneralizedStaticModel('InputName',varargin{1}.InputName,'OutputName',OutNames,'Fun',Fun);
end
MDL = MDL.set('InputMin',min(InMin,[],2),'InputMax',max(InMax,[],2), ...
  	'OutputMin',OutMin,'OutputMax',OutMax,'Parameters',Params,'Free',Free,...
    'NoiseVariance',Nv,'OutputOffset',OutOff);
end