% find_U.m - MMSE receiver
function U = find_U(H,V,sigma2,R,I,K,d)
    U = cell(I,K);
    
    for i = 1:I
        for k = 1:K
            Phi = sigma2 * eye(R);
            
            for j = 1:K
                for l = 1:I
                    Phi = Phi + H{i,k,j} * V{l,j} * V{l,j}' * H{i,k,j}';
                end
            end
            
            % MMSE rcv
            U{i,k} = (Phi \ H{i,k,k} * V{i,k})';
        end
    end
end