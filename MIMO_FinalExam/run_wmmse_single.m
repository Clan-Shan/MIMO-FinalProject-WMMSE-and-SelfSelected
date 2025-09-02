function final_rate = run_wmmse_single(H, P, sigma2, R, I, K, T, d, alpha1, epsilon, max_iter)
    V = cell(I,K);
    for i=1:I
        for k=1:K
            v = randn(T,d)+1i*randn(T,d);
            V{i,k} = sqrt(P/I)*v/norm(v,"fro");
        end
    end
    
    % Run WMMSE algorithm
    rate_old = sum_rate(H,V,sigma2,R,I,K,alpha1);
    iter1 = 1;
    
    while(1)
        U = find_U(H,V,sigma2,R,I,K,d);
        W = find_W(U,H,V,I,K,d);
        V = find_V(alpha1,H,U,W,T,I,K,P);
        rate_new = sum_rate(H,V,sigma2,R,I,K,alpha1);
        iter1 = iter1 + 1;
        if abs(rate_new-rate_old) / rate_old < epsilon || iter1 > max_iter
            break;
        end
        rate_old = rate_new;
    end
    
    final_rate = rate_new;
end