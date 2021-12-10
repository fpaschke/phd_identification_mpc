%1.) Identify simple MISO ARIMAX Model using ls method to estimate initial conditions
rng(0) 
u1 = kron(rand(1,20),ones(1,50))'; 
u2 = kron(rand(1,50),ones(1,20))';
e = .1*randn(size(u1,1),1); u = [u1 u2];
A = conv([1 -.9],[1 -.95]);
C = conv([1 -.7+.4i],[1 -.7-.4i]); B1 = [0 .1]; B2 = [0 .2];
y_ = filter(B1,A,u(:,1)) + filter(B2,A,u(:,2)) + filter(C,A,e);
y = y_ + 0*sin(linspace(0,2*pi,length(u1)))'+0.0;
plot([y y_]);
m1 = idModels.NsfPolyModel(2,1*[1 1],[1 1],2,'Ts',1,'Name','MISO ARIMAX Model');
m1.identify(y,u,'Solver','opti_levmar','CheckGradient','on','PlotZeros',{'C'},'StabilizationMethod','none','Ic','ls',...
    'Hp',1,'Multistep',0,'IntegrateNoise',0,'InitMethod','iv');
assert(all(abs([m1.A.val-A m1.C.val-C m1.B(1,1).val-B1 m1.B(1,2).val-B2])<5e-2));
m1.show('Prediction',y,[u1 u2],'Hp',10)
m1.showBode()
ys = m1.simulate(u);
m1.resample(2)
ys2 = m1.simulate(u(1:2:end,:));
figure(); subplot(2,1,1); plot(y); hold on; plot(ys,'--r'); hold on; plot(1:2:length(u),ys2,'--g'); ylabel('y'); subplot(2,1,2); plot(u); ylabel('u');
close all;

%2.) Estimate AR Model with outputooffset
rng(0)
A = conv([1 -.9],[1 -.95]);
e = randn(1000,1);
y = filter(1,A,e) + 1e3;
m2 = idModels.NsfPolyModel(2,'Name','AR_Model','OutputOffset',5e2);
m2.identify(y,[],'Solver','opti_levmar','CheckGradient','off','EstimateOutputOffset',true);
ys = m2.simulate(zeros(length(y),0));
m2.show('Prediction',y,'Hp',10)
plot([y ys]);
assert(all(abs(m2.A.val-A)<1e-2))
close all;

%3.) Estimate MA Model using stabilization 
rng(0)
C = conv([1 -.9],[1 -.9]);
y = filter(C,1,randn(1000,1));
m3 = idModels.NsfPolyModel(0,[],[],2);
m3.identify(y,[],'Solver','matlab_lsqnonlin','StabilizationMethod','constraint','Ic','backcast');
assert(all(abs(m3.C.val-C)<1e-2))

%4.) Estimate ArArX model 
rng(0) 
u = kron(rand(1,20),ones(1,50))'; 
e = .01*randn(size(u,1),1);
A = conv([1 -.9],1); B = [0 .1]; D = [1 -.8];
y = filter(B,A,u) + filter(1,conv(A,D),e);
m4 = idModels.NsfPolyModel(1,1,1,0,1);
m4.identify(y,u,'Solver','opti_levmar','Ic','backcast');
assert(all(abs([m4.A.val-A m4.D.val-D])<3e-2))

%5.) Estimate overparameterized ArArMaX model
rng(0) 
u = kron(rand(1,200),ones(1,50))'; 
e = .001*randn(size(u,1),1);
A = conv([1 -.9],[1 0]); B = [0 .1]; C = [1 -.7]; D = conv([1 .8],[1]);
y = filter(B,A,u) + filter(C,conv(A,D),e);
m5 = idModels.NsfPolyModel(2,1,1,1,1);
m5.identify(y,u,'Solver','opti_levmar','Hp',1,'Ic','ls');
m5.show('Prediction',y,u,'Hp',10)
assert(all(abs([m5.A.val-A m5.C.val-C m5.D.val-D])<2e-2));


%6.) Identify MISO BJ Model with backcasting of ics
rng(0) 
u1 = kron(rand(1,20),ones(1,50))'; 
u2 = kron(rand(1,50),ones(1,20))';
e = .1*randn(size(u1,1),1); u = [u1 u2];
E = conv([1 -.8],[1 -.85]); C = conv([1 -.5],1); B1 = [0 .1]; B2 = [0 .2]; D = [1 -.9];
y = filter(B1,E,u(:,1)) + filter(B2,E,u(:,2)) + filter(C,D,e);
m6 = idModels.NsfPolyModel(0,[1 1],[1 1],1,1,2);
m6.identify(y,u,'Solver','opti_levmar','ic','backcast','Hp',1,'MultiStep',false,'CheckGradient','on');
dp = [m6.E.val-E m6.C.val-C m6.B(1,1).val-B1 m6.B(1,2).val-B2];
assert(all(abs(dp)<1e-1));

%7.) Identify Full MISO Model with outputoffset
rng(0) 
u1 = kron(rand(1,200),ones(1,50))'; 
u2 = kron(rand(1,500),ones(1,20))';
e = .1*randn(size(u1,1),1); u = [u1 u2];
A = [1 -.9]; E = [1 -.85]; C = [1 -.8]; B1 = [0 .1]; B2 = [0 .2]; D = [1 -.8]; y0 = 10;
y = filter(B1,conv(E,A),u(:,1)) + filter(B2,conv(E,A),u(:,2)) + filter(C,conv(D,A),e) + y0;
m7 = idModels.NsfPolyModel(1,[1 1],[1 1],1,1,1);
m7.identify(y,u,'EstimateOutputOffset',1,'Solver','opti_levmar','CheckGradient','on','InitMethod','iv');
assert(all(abs([y0 - m7.OutputOffset])<1e-1));

% 8.) Identifiy bilinear arx model without noise with outputoffset
clear all; rng(0);
A = conv([1 -.9],1); B = [0 .1]; y0 = 20; p1 = 25;
na = length(A)-1; nk = find(B~=0,1,'first')-1; nb = length(B)-nk; 
F = @(u,p,y) deal(u(:,1).*(p(1) - y),u(:,1));
I(1).fun = F; I(1).parameters = 22; I(1).free = 1; I(1).input_idx = [1]; I(1).output_idx = [1];
u1 = kron(rand(10,1),ones(100,1));
x(1,1) = 0; y = y0;
for k = 2:length(u1)
     [fu,~] = F(u1(k-1:-1:k-nb),p1,y(k-1:-1:k-nb)); 
     x(k,1) = -A(2:end)*x(k-1:-1:k-na) + B(2:end)*fu;
     y(k,1) = y0 + x(k,1);
end
m8 = idModels.NsfPolyModel(1,1,1,'InputNonlinearity',I,'OutputOffset',15);
m8.identify(y,u1,'EstimateOutputOffset',1,'FunTolAbs',1e-20,'Solver','matlab_lsqnonlin'); 
m8.show('Prediction',y,u1,'Hp',10)
figure('color','w'); subplot(2,1,1); plot(y,'-k'); hold on; plot(m8.simulate(u1),'--r'); legend({'y' 'y_s'}); grid on; subplot(2,1,2); plot(u1,'-k'); legend({'u'}); grid on;
assert(all(abs([A - m8.A.val B - m8.B.val p1-m8.InputNonlinearity.parameters{1} y0 - m8.OutputOffset])<1e-3));
close all;

% 9.) Identify Nonlinear ARIMAX Model
rng(0); u = kron(2*rand(20,1),ones(50,1)); N = length(u); e = .01*randn(N,1); 
d = .1*sin(linspace(0,2*pi,N))+0.2; % Some disturbance;
A = conv([1 -.8],[1 -.85]); B = [0 .1]; C = [1 -.85]; 
na = length(A)-1; nk = find(B~=0,1,'first')-1; nb = length(B)-nk; nc = length(C)-1;
alpha = 2;
F = @(u,p,y) deal((u + p(1)*u.^2).*(3-y),... %f
                   u.^2.*(3-y),... %df/dp
                   2*p(1)*u.*(3-y),... %df/du
                   -(u + p(1)*u.^2)); %df/dy
y(1:2,1) = 0; 
for k = 3:N 
    [fu,~,~,~] = F(u(k-1:-1:k-nb),alpha,y(k-1:-1:k-nb)); 
    y(k,1) = -A(2:end)*y(k-1:-1:k-na) + B(2:end)*fu + C*e(k:-1:k-nc) + d(k); 
end
I(1).fun = F; I(1).parameters = [.5]; I(1).free = [true]; I(1).input_idx = 1;
m9 = idModels.NsfPolyModel(2,1,1,1,'InputNonlinearity',I);
m9.identify(y,u,'IntegrateNoise',1,'Solver','matlab_lsqnonlin');
s(1) = subplot(2,1,1); plot(y,'-k'); hold on; plot(m9.simulate(u),'--r'); legend({'y' 'y_{sim}'}); s(2) = subplot(2,1,2); plot(u,'-k'); legend({'u'}); linkaxes(s,'x');
m9.show('Prediction',y,u,'Hp',10)
m9.show('Boxplot',y,u,'Hp',12)
m9.show('Autocovariance',y,u)
Np = 10; ut = linspace(0,1,Np); yt = linspace(0,3,Np);
for k = 1:Np; for l = 1:Np; [f(k,l),~,~,~] = F(ut(k),alpha,yt(l)); [f_hat(k,l),~,~,~] = F(ut(k),m9.InputNonlinearity.parameters{1:end},yt(l)); end; end
figure(); surf(yt,ut,f,'FaceAlpha', 0.2, 'FaceColor', 'b','EdgeColor', 'b'); xlabel('y'); ylabel('u'); zlabel('f(u,y)'); hold on;
surf(yt,ut,f_hat,'FaceAlpha', 0.2, 'FaceColor', 'r','EdgeColor', 'r','LineStyle','--'); xlabel('y'); ylabel('u'); zlabel('f(u,y)'); grid on; legend({'Actual' 'Estimated'},'Location','northeast');
close all;

%10.) Identify ARMA of to low order using k-Step Criterion
close all; clear; 
rng(0); N = 1e3; sigma = 1;
e = sigma*randn(N,1);
A = conv([1 -.85],[1 -.9]);
C = [1 .6];
d = .0*sin(10*pi*linspace(0,1,N)')+0*cos(2*pi*linspace(0,1,N)');
y = filter(C,A,e) + d;
m10 = idModels.NsfPolyModel(1,[],[],1,'OutputName','out');
m10.identify(y,'Hp',1,'solver','opti_levmar','IntegrateNoise',0,'StabilizationMethod','none','initmethod','iv');
figure('color','w'); s = subplot(1,2,1); m10.show('Boxplot',y,'Hp',10,'Axes',s)
m10.identify(y,'Hp',10,'solver','opti_levmar','IntegrateNoise',0,'StabilizationMethod','none','initmethod','iv');
s = subplot(1,2,2); m10.show('Boxplot',y,'Hp',10,'Axes',s)

%11.) Identify simple MISO Hammerstein ARMAX Model 
close all; clear all; rng(0); 
A = conv([1 -.8+.5i],[1 -.8-.5i]); na = length(A)-1; 
u1 = kron(randn(1,20),ones(1,50))'; 
u2 = kron(randn(1,50),ones(1,20))';
u3 = kron(randn(1,10),ones(1,100))';
u = [u1 u2 u3];
B1 = [0 1 1]; B2 = [0 1]; B3 = [0 1]; C = [1 -.8];
e = 1*randn(size(u,1),1);
alpha1 = [1 3]; alpha3 = [1 -2];
I(1).fun = @(u,p) deal(p(1)*u + p(2)*u.^2,[u.^2]); I(1).parameters = [1 2]; I(1).free = [false true]; I(1).input_idx = 1; 
I(2).input_idx = 2;
I(3).fun = @(u,p) deal(p(1)*u + p(2)*u.^2,[u.^2]); I(3).parameters = [1 -1]; I(3).free = [false true]; I(3).input_idx = 3; 
[fu1,~] = I(1).fun(u(:,I(1).input_idx),alpha1); [fu3,~] = I(3).fun(u(:,I(3).input_idx),alpha3); 
y = filter(B1,A,fu1) + filter(B2,A,u(:,2)) + filter(B3,A,fu3) + filter(C,A,e);
subplot(2,1,1); plot(y,'-k'); legend({'y'}); ylabel('y'); subplot(2,1,2); plot(u); ylabel('u');
m11 = idModels.NsfPolyModel(2,[2 1 1],[1 1 1],1,'InputName',{'u1' 'u2' 'u3'},'InputNonlinearity',I);
m11.identify(y,u,'Solver','matlab_lsqnonlin','CheckGradient','on');
m11.show('Prediction',y,u,'Hp',10);
figure(); ut = repmat((0:.01:3)',1,3);
[fu1_r,~] = I(1).fun(ut,alpha1); [fu3_r,~] = I(3).fun(ut,alpha3); plot(ut(:,1),[fu1_r fu3_r]); hold on; grid on;
[fu1_t,~] = I(1).fun(ut,cell2mat(m11.InputNonlinearity(1).parameters)); [fu3_t,~] = I(3).fun(ut,cell2mat(m11.InputNonlinearity(3).parameters)); plot(ut(:,1),[fu1_t fu3_t],'--'); legend({'Actual f1' 'Actual f3' 'Estimated f1' 'Estimated f3'});
close all;

% 12.) Ident 2 Models merge and simulate them
rng(0);
u1 = kron(rand(1,20),ones(1,50))'; 
u2 = kron(rand(1,50),ones(1,20))';
e = .1*randn(size(u1,1),1); u = [u1 u2];
A1 = conv([1 -.9],[1 -.95]);
C1 = conv(1,1); B11 = [0 .1]; B12 = [0 .2];
y1 = filter(B11,A1,u(:,1)) + filter(B12,A1,u(:,2)) + filter(C1,A1,e);
m12_1 = idModels.NsfPolyModel(2,[1 1],[1 1],2,'InputName',{'u1' 'u2'},'OutputName','y1');
m12_1.identify(y1,u,'Solver','opti_levmar');
m12_1.show('Boxplot',y1,u,'Hp',20)
A2 = conv([1 -.95],1);
C2 = 1; B21 = [0 .1]; B22 = [0 .2];
y2 = filter(B21,A2,u(:,1)) + filter(B22,A2,u(:,2)) + filter(C2,A2,e);
m12_2 = idModels.NsfPolyModel(1,[1 1],[1 1],0,'InputName',{'u3' 'u4'},'OutputName','y2');
m12_2.identify(y2,u,'Solver','opti_levmar');
m12_2.show('Boxplot',y2,u,'Hp',20)
m12 = m12_1.merge(m12_2);
m12.show('Boxplot',[y1 y2],[u u],'Hp',20);

% 13.) Ident 2 Models merge and simulate them
rng(0); Ts = 1;
u1 = kron(rand(1,20),ones(1,50))'; 
u2 = kron(rand(1,50),ones(1,20))';
e = .1*randn(size(u1,1),1); u = [u1 u2];
A1 = conv([1 -.9],1);
C1 = [1 .8]; B11 = [0 .1]; B12 = [0 .2];
y1 = filter(B11,A1,u(:,1)) + filter(B12,A1,u(:,2)) + filter(C1,A1,e);
m13_1 = idModels.NsfPolyModel(1,[1 1],[1 1],1,'InputName',{'u1' 'u2'},'OutputName','y1');
m13_1.identify(y1,u,'Solver','opti_levmar');
m13_1.show('Boxplot',y1,u,'Hp',20)
% m13_1.resample(Ts);
% m13_1.show('Boxplot',y1(1:Ts:end),u(1:Ts:end,:),'Hp',10)
A2 = conv([1 -.95],1);
C2 = 1; B21 = [0 .1]; B22 = [0 .2];
y2 = filter(B21,A2,u(:,1)) + filter(B22,A2,u(:,2)) + filter(C2,A2,e);
m13_2 = idModels.NsfPolyModel(1,[1 1],[1 1],0,'InputName',{'u3' 'u4'},'OutputName','y2');
m13_2.identify(y2,u,'Solver','opti_levmar');
m13_2.show('Boxplot',y2,u,'Hp',20)
% m13_2.resample(Ts);
% m13_2.show('Boxplot',y2(1:Ts:end),u(1:Ts:end,:),'Hp',10)
m13 = m13_1.merge(m13_2);
m13.show('Boxplot',[y1(1:Ts:end) y2(1:Ts:end)],[u(1:Ts:end,:) u(1:Ts:end,:)],'Hp',10);

%14.) Identify MIMO ARMAX Model with outputoffset
A = {[-.9 .5; 0 -.8]}; ny = size(A{1},2); B = {zeros(ny,2) eye(ny,2)}; nu = size(B{1},2); C = {[.8 .0; .0 .7]};
rng(0); Np = 20; Ns = 50; std_e = [.1 .05]; u = kron(rand(Np,nu),ones(Ns,1)); e = std_e.*randn(Np*Ns,ny);
y_ = idModels.alg.lti.lsim_poly(A,B,[],u);
y = idModels.alg.lti.lsim_poly(A,B,C,u,e) + [5 15];
subplot(2,1,1); plot(y(:,1),'-k'); hold on; plot(y_(:,1),'--k'); yyaxis right;
plot(y(:,2),'-b'); hold on; plot(y_(:,2),'--b'); subplot(2,1,2); plot(u);
m14 = idModels.NsfPolyModel([1 1; 0 1],[1 1; 1 1],[1 1; 1 1],[1 1],'OutputOffset',[5 10]);
m14.identify(mat2cell(y,[500 500],2),mat2cell(u,[500 500],2),'EstimateOutputOffset',[0 1],'Solver','opti_levmar','CheckGradient','on','MultiStep',false,'Hp',1*[1 1],'Ic','ls');
m14.show('Boxplot',y,u,'Hp',10);
m14.show('Prediction',y,u,'Hp',50);
m14.showBode()

%15.) Identify MIMO BJ/ARMAX Model with multistep criterion
E1 = [1 -.9]; E2 = [1 -.85]; B11 = [0 1]; B12 = [0 1]; B21 = [0 1]; B22 = [0 1]; C1 = [1 .7]; C2 = [1 .7];
rng(0); Np = 20; Ns = 50; std_e = [.1 .1]; u = kron(rand(Np,2),ones(Ns,1)); e = std_e.*randn(Np*Ns,2);
y1 = filter(B11,E1,u(:,1)) + filter(B12,E1,u(:,2)) + filter(C1,E1,e(:,1));
y2 = filter(B21,E2,u(:,1)) + filter(B22,E2,u(:,2)) + filter(C2,E2,e(:,2));
subplot(2,1,1); plot(y1,'-k'); yyaxis right; plot(y2,'-b'); subplot(2,1,2); plot(u);
m15 = idModels.NsfPolyModel(zeros(2),[1 1; 1 1],[1 1; 1 1],[1 1],[1 1],[1 1]);
m15.identify([y1 y2],u,'Solver','matlab_lsqnonlin','CheckGradient','off','MultiStep',false,'Hp',[3 4],'Ic','ls');
m15.show('Boxplot',[y1 y2],u,'Hp',10);
m15.show('Prediction',[y1 y2],u,'Hp',10);
m15 = idModels.NsfPolyModel(eye(2),[1 1; 1 1],[1 1; 1 1],[1 1]);
m15.identify([y1 y2],u,'Solver','opti_levmar','CheckGradient','on','MultiStep',false,'Hp',[4 6],'Ic','ls');
m15.show('Boxplot',[y1 y2],u,'Hp',10);
m15.show('Prediction',[y1 y2],u,'Hp',10);

%16.) Identify ARMAX with Prefilter
rng(0) 
u = kron(rand(1,20),ones(1,50))'; 
e = .1*randn(size(u,1),1);
A = conv(1,[1 -.9]);
B = [0 1];
C = conv(1,[1 .8]); 
Zp = [1 -1];
Np = [1 .7];
y = filter(B,A,u) + filter(conv(C,Np),conv(A,Zp),e);
m16 = idModels.NsfPolyModel(1,1,1,1,'Ts',1,'Name','SISO ARIMAX Model','PreFilter',{Zp Np});
m16.identify(y,u,'Hp',10,'Solver','opti_levmar','CheckGradient','on');
assert(all(abs([m16.A.val-A m16.B.val-B m16.C.val-C])<5e-2));
m16.show('Prediction',y,u,'Hp',10)

%17.) Identify BJ with Prefilter
rng(0) 
u = kron(rand(1,200),ones(1,50))'; 
e = .1*randn(size(u,1),1);
B = [0 1];
C = conv(1,[1 .8]); 
D = conv(1,[1 -.8]); 
E = conv(1,[1 -.9]);
Zp = [1 -1];
Np = [1 .7];
y = filter(B,E,u) + filter(conv(C,Np),conv(D,Zp),e);
m17 = idModels.NsfPolyModel(0,1,1,1,1,1,'Ts',1,'Name','SISO BJ','PreFilter',{Zp Np});
m17.identify(y,u,'Hp',5,'Solver','opti_levmar','CheckGradient','on');
assert(all(abs([m17.B.val-B m17.C.val-C m17.D.val-D m17.E.val-E])<5e-2));
m17.show('Prediction',y,u,'Hp',10);

%18. Identify ARMA Model with factorized A and outputoffset using multistep criterion
A = conv([1 -.9],[1 -.8]);
C = conv(1,[1 .8]);
rng(0); e = randn(1000,1);
y = filter(C,A,e) + 100;
m18 = idModels.NsfPolyModel(2,[],[],1);
m18.factorize('A');
m18.identify(y,'Solver','matlab_lsqnonlin','CheckGradient','on','EstimateOutputOffset',true,'Hp',4,'Multistep',true);
assert(all([A - idModels.util.defactorizePoly(m18.A.val) C - m18.C.val]<3e-2));

%19. Identify ARMAX Model with factorized A and B Prefilter
A = conv([1 -.7],[1 -.8]);
B = [0 0 2*conv([1 -.9],[1 -.5])];
C = conv(1,[1 .8]);
Zp = [1 -.5];
Np = [1 -.0];
rng(0); e = randn(10000,1); u = rand(length(e),1);
y = filter(B,A,u) + filter(C,conv(A,Zp),e);
m19 = idModels.NsfPolyModel(2,3,2,1,'PreFilter',{Zp Np});
m19.factorize('A'); m19.A.val = [1 .5 .4];
m19.factorize('B'); m19.B.val = [1 Inf Inf .4 .6];
m19.factorize({'A' 'B'});
m19.identify(y,u,'Solver','opti_levmar','CheckGradient','on');
m19.showBode
assert(all([A - idModels.util.defactorizePoly(m19.A.val) B - idModels.util.defactorizePoly(m19.B.val) C - m19.C.val]<5e-2));

%20.) Identify MIMO ARMAX Model with factorized A and Outputoffset
A = {[-.9 .5; 0 -.8]}; ny = size(A{1},2); B = {zeros(ny,2) eye(ny,2)}; nu = size(B{1},2); C = {[.8 .0; .0 .7]};
rng(0); Np = 20; Ns = 50; std_e = [.1 .05]; u = kron(rand(Np,nu),ones(Ns,1)); e = std_e.*randn(Np*Ns,ny);
y_ = idModels.alg.lti.lsim_poly(A,B,[],u);
y = idModels.alg.lti.lsim_poly(A,B,C,u,e) + [5 15];
subplot(2,1,1); plot(y(:,1),'-k'); hold on; plot(y_(:,1),'--k'); yyaxis right;
plot(y(:,2),'-b'); hold on; plot(y_(:,2),'--b'); subplot(2,1,2); plot(u);
m20 = idModels.NsfPolyModel([1 1; 1 1],[1 1; 1 1],[1 1; 1 1],[1 1],'OutputOffset',[5 10]);
m20.factorize('A');
m20.A = {[1 .6] [.1 Inf]; [.1 Inf] [1 .5]};
m20.identify(mat2cell(y,[500 500],2),mat2cell(u,[500 500],2),'EstimateOutputOffset',[0 1],'Solver','opti_levmar','CheckGradient','on');
m20.show('Boxplot',y,u,'Hp',10);
m20.show('Prediction',y,u,'Hp',10);

%21. Identify factorized BJ Model with Prefilter
E = conv([1 -.9],1);
D = [1 -.7];
B = [0 0 2*conv([1 .6],[1 -.4])];
C = conv(1,[1 .8]);
Zp = [1 -1];
Np = [1 -.0];
rng(0); e = randn(1e5,1); u = rand(length(e),1);
y = filter(B,E,u) + filter(C,conv(D,Zp),e);
m21 = idModels.NsfPolyModel(0,3,2,1,1,1,'PreFilter',{Zp Np});
m21.factorize();
m21.identify(y,u,'Solver','matlab_lsqnonlin','CheckGradient','off','Ic','ls','Hp',1);
m21.showBode('ShowH',true)
m21.defactorize();
assert(all([E - m21.E.val B - m21.B.val C - m21.C.val D - m21.D.val]<5e-2));

%22. Identify factorized BJ Model
B = [0 0 2*conv([1 .5],[1 -.5])];
C = conv([1 -.5],[1 .8]);
D = conv([1 -.8],[1 .5]);
E = [1 -.9];
rng(0); e = randn(10000,1); u = rand(length(e),1);
y = filter(B,E,u) + filter(C,D,e);
m22 = idModels.NsfPolyModel(0,3,2,2,2,1); 
m22.factorize();
m22.B.val = [1 Inf Inf .6 -.6];
m22.C.val = [1 .6 -.7];
m22.D.val = [1 .7 -.6];
m22.E.val = [1 .8];
m22.identify(y,u,'Solver','opti_levmar','CheckGradient','on','Ic','backcast','Hp',1);
m22.showBode('ShowH',true)
assert(all(abs([B - idModels.util.defactorizePoly(m22.B.val) C - idModels.util.defactorizePoly(m22.C.val) D - idModels.util.defactorizePoly(m22.D.val)])<6e-2));

%23. Identify Hammerstein Model with shared parameter
clear; rng(0); u = rand(1e3,1);
A = [1 -.9]; B1 = [0 1]; B2 = [0 0 1];
p1 = 1; p2 = 2;  p3 = 3;
I(1).fun = @(u,p) deal(1 + p(1)*u.^2 + p(2)*u.^3,[u.^2 u.^3]);
I(1).input_idx = 1;
I(1).parameters = {'2_2' p1+1}; 
I(1).free = 1;
I(2).fun = @(u,p) deal(1 + p(1)*log(u) + p(2)*u.^2,[log(u) u.^2]); 
I(2).input_idx = 1;
I(2).parameters = {p2 p3+1}; 
I(2).free = [0 1];
[f1 , ~] = I(1).fun(u,[p3 p1]); % eval f1
[f2 , ~] = I(2).fun(u,[p2 p3]); % eval f2
y = filter(B1,A,f1) + filter(B2,A,f2); % create output;
m23 = idModels.NsfPolyModel(1,[1 2],[1 1],'InputNonlinearity',I); 
m23.identify(y,u,'solver','matlab_lsqnonlin','CheckGradient','on');
assert(all(abs([m23.InputNonlinearity(1).parameters{2}-p1 m23.InputNonlinearity(2).parameters{1}-p2 m23.InputNonlinearity(2).parameters{2}-p3])<1e-6));

%24. Identify MIMO ARMAX Model with Prefilter
A = {[-.9 .5; 0 -.8]}; ny = size(A{1},2); B = {zeros(ny,2) eye(ny,2)}; nu = size(B{1},2); C = {[.8 .0; .0 .7]};
rng(0); Np = 20; Ns = 50; std_e = [.1 .05]; u = kron(rand(Np,nu),ones(Ns,1)); e = 1*std_e.*randn(Np*Ns,ny);
y_ = idModels.alg.lti.lsim_poly(A,B,[],u);
y = idModels.alg.lti.lsim_poly(A,B,C,u,e) + [5 15];
subplot(2,1,1); plot(y(:,1),'-k'); hold on; plot(y_(:,1),'--k'); yyaxis right;
plot(y(:,2),'-b'); hold on; plot(y_(:,2),'--b'); subplot(2,1,2); plot(u);
Pf.num = [1 -1];
Pf.den = [1 -.9];
m24 = idModels.NsfPolyModel([1 1; 0 1],[1 1; 1 1],[1 1; 1 1],[1 1],'OutputOffset',[5 10],'PreFilter',Pf);
m24.identify(y,u);
m24.show('Prediction',y,u,'Hp',50);
m24.showBode()
