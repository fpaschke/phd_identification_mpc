function [G_num,G_den,H_num,H_den] = calcGH(obj,varargin)
%CALCGH Returns deterministic Part G and stochastic part H of lti part of the model, such that
% y = G(q^-1) g(u(t),y(t)) + H(q^-1) e(t).
%
% [G_num,G_den,H_num,G_den] = calcGH(obj,varargin)
% obj [NsfPolyModel]: The Model object.
% G_num [ny x nv cell of double vectors]: Numerator part of G.
% G_den [ny x nv cell of double vectors]: Denumerator part of G.
% H_num [ny x 1 cell of double vectors]: Numerator part of H.
% H_num [ny x 1 cell of double vectors]: Denumerator part of H.
% varargin: Optinal Name Valiue Pairs
%   AbsorbNoiseVariance [logical]: Absorb NoiseVariance to Model such that e has Variance = 1 [Def. true]

p = inputParser();
addParameter(p,'AbsorbNoiseVariance',true,@(x) x == 1 || x == 0);
addParameter(p,'TolCancelPolesZeros',[],@(x) isempty(x) || isscalar(x) && x>=0 && isreal(x));
addParameter(p,'TolDiscardPolesZeros',[],@(x) isempty(x) || isscalar(x) && x>=0 && isreal(x));
addParameter(p,'ReduceGonly',[],@(x) isempty(x) || isscalar(x) && x>=0 && isreal(x));
parse(p,varargin{:});

[A,B,C,D,E,F] = idModels.util.polyStruct2cell(obj.A,obj.B,obj.C,obj.D,obj.E,obj.F,'DefactorizePoly',true,'ImagTol',obj.ImagTol);
isDiagA = isdiag(double(cellfun(@(x) ~all(x==0) && ~isempty(x),A)));
assert(obj.IsParameterized && all(~isnan(obj.NoiseVariance(:))),'The Model is not fully parameterized i. e. it contains NaN parameters or any value of the noise covariance (NoiseVariance) is NaN!');
nu = size(obj.Nb,2); % #ins to linear block
ny = obj.OutputDimension;
G_num = cell(obj.OutputDimension,size(obj.B,2)); G_den = G_num;
H_num = cell(obj.OutputDimension,obj.OutputDimension); H_den = H_num;

if isDiagA
    for no = 1:ny
        % Deterministic part
        for ni = 1:nu
            G_num{no,ni} = B{no,ni};
            G_den{no,ni} = conv(conv(A{no,no},E{no}),F{no,ni});
        end
        % Noise
        if nargout > 2
            for no2 = 1:ny
                if no == no2
                    H_num{no,no} = conv(C{no},obj.PreFilter(no).den);
                    H_den{no,no} = conv(conv(A{no,no},D{no}),obj.PreFilter(no).num);
                    if p.Results.AbsorbNoiseVariance
                       H_num{no,no} = H_num{no,no}*sqrt(obj.NoiseVariance(no,no));
                    end
                else
                    H_num{no,no2} = 0;
                    H_den{no,no2} = 1;
                end
            end
        end
    end
else
    for no=1:ny 
        C{no} = conv(C{no},obj.PreFilter(no).den);
        D{no} = conv(D{no},obj.PreFilter(no).num);
        F(no,:) = cellfun(@(fi) conv(fi,E{no}),F(no,:),'UniformOutput',false);
    end
    GH = tf(idpoly(A,B,C,D,F,obj.NoiseVariance,obj.Ts,'TimeUnit',obj.TimeUnit),'augmented');
    if isfield(GH.InputGroup,'Measured')
     	G_num = GH.Numerator(:,GH.InputGroup.Measured);
        G_den = GH.Denominator(:,GH.InputGroup.Measured);
    else
        G_num = [];
        G_den = [];
    end
   	H_num = GH.Numerator(:,GH.InputGroup.Noise);
    H_den = GH.Denominator(:,GH.InputGroup.Noise);
end
end


