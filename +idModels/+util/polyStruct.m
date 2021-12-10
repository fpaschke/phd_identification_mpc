function P = polyStruct(n,type,nd)
% Creates size(n)-Polystructure P with fields val and free which represent a matrix polynomial 
%
%  P(q^-1) = P_0 + P_1q^-1 + ... 
% 
% P = polyStruct(n,type,nd)
% P [struct]: Polystructure
% n [pos. int matrix]: defines orders of polynomials
% nd [pos. int matrix]: defines delays of polynomials
% type ['a' 'b' or 'f']:    'a': P_0(0) will be  eye(size(n)) -> use for init of A,C,D,E polynomials
%                           'b': will set nd(r,c) delays tp P -> use for init of B polynomials

P = repmat(struct('val',[],'free',[],'factorized',false),size(n,1),size(n,2));
if ~isstruct(n) % init Structures
    for r = 1:size(n,1)
        for c = 1:size(n,2)
            P(r,c).factorized = false;
            switch lower(type)
                case 'a'
                    if r == c
                        P(r,c).val = [1 NaN(1,n(r,c))];
                    else
                        P(r,c).val = [0 NaN(1,n(r,c))];
                    end
                    P(r,c).free = [false true(1,n(r,c))]; 
                case 'b'
                    P(r,c).free = [false(1,nd(r,c)) true(1,n(r,c) - nd(r,c) + 1)]; 
                    P(r,c).val =  [zeros(1,nd(r,c)) NaN(1,n(r,c) - nd(r,c) + 1)];
                case 'f'
                    if r == c
                        P(r,c).val = [1 NaN(1,n(r,c))];
                        P(r,c).free = [false true(1,n(r,c))]; 
                    else
                        P(r,c).val = 0;
                        P(r,c).free = false; 
                    end
            end
        end
    end
else
    P = n;
end
if diff(size(P))~=0 && ~strcmpi(type,'b') % diagonalize
    for r = 1:length(P)
       for c = 1:length(P)
           if r == c
               X(r,c) = P(r);
           else
               X(r,c).val(1) = 0;
               X(r,c).free(1) = false;
               X(r,c).factorized = false;
           end
       end
    end
    P = X;
end
end

