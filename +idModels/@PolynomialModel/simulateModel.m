function y = simulateModel(obj,u)
% SIMULATEMODEL Simulates polynomial Model for given Inputs u.
%
% y = simulateModel(obj,u)
% obj [PolynomialModel]: Model Object
% u [N x nu double]: Inputs
% y [N x ny double]: Simulated Output

%assert(obj.OutputDimension == 1);

y = NaN(size(u,1),obj.OutputDimension);
PHI = getRegressorMatrices(obj,u);
for no = 1:obj.OutputDimension
    y(:,no) = cell2mat(PHI(no,:))*[obj.Coefficients(no,:).val]' + obj.OutputOffset(no);
end
end