% find_W.m - Weight matrix update
function W = find_W(U,H,V,I,K,d)
    W = cell(I,K);
    
    for i = 1:I
        for k = 1:K
            E = eye(d) - U{i,k} * H{i,k,k} * V{i,k};

            try
                W{i,k} = inv(E);
            catch
                % E is singular, use pseudo-inverse
                W{i,k} = pinv(E);
            end
        end
    end
end