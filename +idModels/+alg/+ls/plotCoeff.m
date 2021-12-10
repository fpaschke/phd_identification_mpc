function plotCoeff(varargin)
%PLOTCOEFF Summary of this function goes here
%   Detailed explanation goes here
    f = figure('Color','w'); 
    tabgp = uitabgroup(f);
    names = {'A' 'B' 'C' 'D' 'E' 'F'};
    for np = 1:length(varargin)
        if ~isempty(varargin{np})
            tab1 = uitab(tabgp,'Title',names{np});
            nsp = 1;
            [nr,nc] = size(varargin{np}{1});
            for r = 1:nr
                for c = 1:nc
                    Pi = cell2mat(cellfun(@(ai) ai{r,c},varargin{np},'UniformOutput',false));
                    subplot(nr,nc,nsp,'Parent',tab1); 
                    if np == 2 % only B isnt Monic
                        plot(Pi(:,1:end)); legend(arrayfun(@(i) [lower(names{np}) '_' num2str(i)],0:size(Pi,2)-1,'UniformOutput',false));
                    else
                        plot(Pi(:,2:end)); legend(arrayfun(@(i) [lower(names{np}) '_' num2str(i)],1:size(Pi,2)-1,'UniformOutput',false));
                    end
                    grid on;
                    nsp = nsp + 1; 
                end
            end
        end
    end
end

