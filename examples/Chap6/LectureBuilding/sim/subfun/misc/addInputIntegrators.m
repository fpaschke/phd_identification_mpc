function M = addInputIntegrators(mdl,ins)
% Adds integrators to specified input channels.

ni = length(ins);
ny = mdl.OutputDimension;
f = mdl.InputNonlinearity;
nv = length(f);
P0 = struct('val',0,'std',NaN,'free',false,'factorized',false); % P = 0
P1 = struct('val',1,'std',NaN,'free',false,'factorized',false); % P = 1
Pi = struct('val',[1 -1],'std',[NaN NaN],'free',false(1,2),'factorized',false); % Integrator
A = [ mdl.A  repmat(P0,ny,ni);
      repmat(P0,ni,ny) repmat(P0,ni,ni)];
B = [ mdl.B repmat(P0,ny,ni);
      repmat(P0,ni,nv) repmat(P0,ni,ni)];
C = [ mdl.C; repmat(P1,ni,1)];  
D = [ mdl.D; repmat(P1,ni,1)];  
E = [ mdl.E; repmat(P1,ni,1)];  
F = [ mdl.F repmat(P1,ny,ni);
      repmat(P1,ni,nv) repmat(P1,ni,ni)]; 
Out = [mdl.OutputName; mdl.InputName(ins)];
In = mdl.InputName;
In(ins) = cellfun(@(n) ['Delta ' n],In(ins),'UniformOutput',0);
Lambda = [mdl.NoiseVariance zeros(ny,ni);
          zeros(ni,ny) zeros(ni,ni)];
Pf = [mdl.PreFilter;
      repmat(struct('num',1,'den',1),ni,1)];

for nf = 1:nv % replace input_idx by output_idx 
    for i = ins 
        ix = f(nf).input_idx == i;
        f(nf).input_idx(ix) = [];
        if any(ix)
            f(nf).output_idx = [f(nf).output_idx ny+find(i==ins)];
        end
    end
end

for i = 1:ni
    A(ny+i,ny+i) = P1; % Set 1 to diagonals
    B(ny+i,nv+i) = P1; % Set 1
    F(ny+i,nv+i) = Pi; % Set Integrators
    f = [f; struct('fun',[],'input_idx',ins(i),'parameters',[],'free',[],'output_idx',[],'A',[],'b',[],'std',[])];
end  
M = idModels.NsfPolyModel(zeros(size(A)),zeros(size(B)),zeros(size(B)),zeros(size(C)),zeros(size(D)),zeros(size(E)),zeros(size(F)),...
        'A',A,'B',B,'C',C,'D',D,'E',E,'F',F,'NoiseVariance',Lambda,'InputName',In,'OutputName',Out,...
    	'Ts',mdl.Ts,'InputNonlinearity',f,'PreFilter',Pf,'TimeUnit',mdl.TimeUnit,...
        'OutputOffset',[mdl.OutputOffset; zeros(ni,1)],'Name',mdl.Name); 

end

