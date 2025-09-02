% find_V.m - Transmit beamformer update
function V = find_V(alpha1,H,U,W,T,I,K,P)
    V = cell(I,K);
    
    for k = 1:K
        A_k = zeros(T,T);
        B_k = zeros(T,I*size(W{1,k},1));
        
        col_idx = 1;
        for i = 1:I
            A_k = A_k + alpha1(i,k) * H{i,k,k}' * U{i,k}' * W{i,k} * U{i,k} * H{i,k,k};
            d_size = size(W{i,k},1);
            B_k(:,col_idx:col_idx+d_size-1) = alpha1(i,k) * H{i,k,k}' * U{i,k}' * W{i,k};
            col_idx = col_idx + d_size;
        end
        
        for i = 1:I
            for j = 1:K
                if j ~= k
                    A_k = A_k + alpha1(i,j) * H{i,j,k}' * U{i,j}' * W{i,j} * U{i,j} * H{i,j,k};
                end
            end
        end
        
        for i = 1:I
            % Use water-filling or simple scaling to satisfy power constraint
            if rank(A_k) == size(A_k,1)
                d_size = size(W{i,k},1);
                start_idx = (i-1)*d_size + 1;
                end_idx = i*d_size;
                
                V_temp = A_k \ B_k(:,start_idx:end_idx);
                
                % Scale to satisfy power constraint
                current_power = norm(V_temp,'fro')^2;
                if current_power > P/I
                    V{i,k} = sqrt(P/I) * V_temp / sqrt(current_power);
                else
                    V{i,k} = V_temp;
                end
            else
                % If A_k is singular
                v = randn(T,size(W{i,k},1))+1i*randn(T,size(W{i,k},1));
                V{i,k} = sqrt(P/I) * v / norm(v,"fro");
            end
        end
    end
end