% sum_rate.m
function system_rate = sum_rate(H,V,sigma2,R,I,K,alpha1)
    system_rate = 0;
    
    for i = 1:I
        for k = 1:K
            Phi = sigma2 * eye(R);
            
            for j = 1:K
                for l = 1:I
                    Phi = Phi + H{i,k,j} * V{l,j} * V{l,j}' * H{i,k,j}';
                end
            end
            
            % Calculate rate for user i_k
            numerator = det(Phi);
            denominator = det(Phi - H{i,k,k} * V{i,k} * V{i,k}' * H{i,k,k}');
            
            if denominator > 0
                rate_ik = real(log2(numerator / denominator));
            else
                rate_ik = 0;
            end
            
            system_rate = system_rate + alpha1(i,k) * rate_ik;
        end
    end
end
