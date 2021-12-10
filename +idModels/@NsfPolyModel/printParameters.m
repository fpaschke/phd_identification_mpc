function printParameters(obj,varargin)
%PRINTPARAMETERS Prints estimated Parameters and their Standarddeviation.
%Model needs to be parameterized.

p = inputParser();
addParameter(p,'PrintFixedZeros',false,@(x) x == 0 || x ==1);
addParameter(p,'Ndigits',4,@(x) mod(x,1)==0 && x>0);
addParameter(p,'PrintRelativeDeviation',true,@(x) x == 0 || x ==1);
parse(p,varargin{:});

assert(obj.IsParameterized,'Model is not parameterized!')
locPrintPoly(obj,'A',p.Results);
locPrintPoly(obj,'B',p.Results);
locPrintPoly(obj,'C',p.Results);
locPrintPoly(obj,'D',p.Results);
locPrintPoly(obj,'E',p.Results);
locPrintPoly(obj,'F',p.Results);
locPrintInputNonlinearity(obj,p.Results)
end

function locPrintInputNonlinearity(obj,opt)
    for ni = 1:length(obj.InputNonlinearity)
        if ~isempty(obj.InputNonlinearity(ni).fun)
            if ~isempty(obj.InputNonlinearity(ni).output_idx)
                str = 'INPUTNONLINEARITY WITH OUTPUT FEEDBACK'; 
                if ~isempty(obj.InputNonlinearity(ni).input_idx)
                    str2 = 'u,y,eta';
                else
                    str2 = 'u,eta';
                end
            else
                str = 'HAMMERSTEIN INPUTNONLINEARITY'; str2 = 'u,eta';
            end
            fprintf('%s f_%i(%s_%i): \n',str,ni,str2,ni);
            nfp = 1;
            for np = 1:length(obj.InputNonlinearity(ni).parameters)
                if ~ischar(obj.InputNonlinearity(ni).parameters{np})
                    if obj.InputNonlinearity(ni).free(nfp)
                        str = ['+- ' num2str(100*obj.InputNonlinearity(ni).std(nfp)/abs(obj.InputNonlinearity(ni).parameters{np}),opt.Ndigits) '%']; 
                    else
                        str = '(fixed)';
                    end
                    fprintf('eta_%i,%i = %s %s \n',ni,np,num2str(obj.InputNonlinearity(ni).parameters{np},opt.Ndigits),str);
                    nfp = nfp + 1;
                end
            end
        end
    end
end

function locPrintPoly(obj,P,opt)
    Poly = obj.(P);
    [Nr,Nc] = size(Poly);
    if ~strcmpi(P,'B') && all(all(obj.(['N' lower(P)])==0)) || strcmpi(P,'B') && all(all(obj.(['N' lower(P)])==-1)) % Do not print 
        return;
    end
    for r = 1:Nr
        for c = 1:Nc
            if isequal(Poly(r,c).val,0)
                continue;
            end
            if Poly(r,c).factorized
                str = 'FACTORIZED POLYNOMIAL '; 
            else
                str = 'POLYNOMIAL ';
            end
            fprintf([str '%s(%i,%i): \n'],P,r,c);
            for np = 1:length(Poly(r,c).val)
                if Poly(r,c).factorized
                    if np == 1
                        p = 'K = ';
                    else
                        p = [lower(P) '_' num2str(np-1) ',0 = '];
                    end
                else
                    if opt.PrintFixedZeros==false && Poly(r,c).val(np)==0
                        continue;
                    end
                    p = [lower(P) '_' num2str(np-1) ' = '];
                end
                p = [p num2str(Poly(r,c).val(np),opt.Ndigits)];
                if Poly(r,c).free(np)
                    if opt.PrintRelativeDeviation
                        p = [p ' +- ' num2str(Poly(r,c).std(np)/abs(Poly(r,c).val(np))*100,opt.Ndigits) '%%'];
                    else
                        p = [p ' +- ' num2str(Poly(r,c).std(np),opt.Ndigits)];
                    end
                else
                    if Poly(r,c).factorized && np ~= 1 && isinf(Poly(r,c).val(np))
                        p = [p ' (fixed delay)'];
                    else
                        p = [p ' (fixed)'];
                    end
                end
                fprintf([p '\n']);
            end
        end
    end
end