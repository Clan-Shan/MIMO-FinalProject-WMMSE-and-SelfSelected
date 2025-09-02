clear;
clc;
close all;

K = 4; % BS
T = 1; % transmit antennas
R = 1; % receive antennas
epsilon = 1e-3;
sigma2 = 1; % noise power
I = 4; % users per BS
alpha1 = ones(I,K); % weight
d = 1; % data stream
snr = 30;
P = db2pow(snr)*sigma2; % transimt power
H = cell(I,K,K); % channel
max_iter = 100;

for i=1:I
    for k = 1:K
        for j=1:K
           H{i,k,j}=sqrt(1/2)*(randn(R,T)+1i*randn(R,T));
        end
    end
end
rate = []; 

% U = cell(I,K);
% U(:)={zeros(R,d)};

V = cell(I,K); 
for i=1:I
    for k=1:K
        v = randn(T,d)+1i*randn(T,d);
        V{i,k}=sqrt(P/I)*v/norm(v,"fro");
    end
end 

rate_old = sum_rate(H,V,sigma2,R,I,K,alpha1);
rate = [rate rate_old];

iter1 = 1;
while(1)
    U = find_U(H,V,sigma2,R,I,K,d); 
    W = find_W(U,H,V,I,K,d); 
    V = find_V(alpha1,H,U,W,T,I,K,P); 
    rate_new = sum_rate(H,V,sigma2,R,I,K,alpha1);
    rate = [rate rate_new];
    iter1 = iter1 + 1;
    if abs(rate_new-rate_old) / rate_old < epsilon || iter1 > max_iter
        break;
    end
    rate_old = rate_new;
end

plot(0:iter1-1,rate,'r-o')
grid on
xlabel('Iterations')
ylabel('Sum rate (bits per channel use)')
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1)
% title('MIMO-IFC, SNR=30, K=4, T=4, R=2, \epsilon=1e-3','Interpreter','tex')
title('SISO-IFC, SNR=30, K=4, T=1, R=1, \epsilon=1e-3','Interpreter','tex')


run_snr_analysis(3, 2, 2, 2, 2, 'MIMO');
run_snr_analysis(3, 1, 1, 2, 1, 'SISO');

function run_snr_analysis(K, T, R, I, d, case_name)
    epsilon = 1e-3;
    sigma2 = 1;
    max_iter = 100;
    alpha1 = ones(I,K);
    snr_range = 5:5:30;
    num_realizations = 100;

    rate_snr = zeros(size(snr_range));
    rate_snr_int10 = zeros(size(snr_range));

    for idx = 1:length(snr_range)
        snr_db = snr_range(idx);
        P = db2pow(snr_db) * sigma2;

        rates_single = zeros(1, num_realizations);
        rates_multi = zeros(1, num_realizations);

        fprintf('Processing SNR = %d dB (%s)\n', snr_db, case_name);

        parfor real = 1:num_realizations
            % Generate random channel
            H = cell(I,K,K);
            for i=1:I
                for k = 1:K
                    for j=1:K
                        H{i,k,j} = sqrt(1/2)*(randn(R,T)+1i*randn(R,T));
                    end
                end
            end

            % Single run
            rates_single(real) = run_wmmse_single(H, P, sigma2, R, I, K, T, d, alpha1, epsilon, max_iter);

            % 10 random initializations
            best_rate = 0;
            for init = 1:10
                rate = run_wmmse_single(H, P, sigma2, R, I, K, T, d, alpha1, epsilon, max_iter);
                if rate > best_rate
                    best_rate = rate;
                end
            end
            rates_multi(real) = best_rate;
        end

        rate_snr(idx) = mean(rates_single);
        rate_snr_int10(idx) = mean(rates_multi);
    end

    % Plot results
    figure;
    plot(snr_range,rate_snr_int10,'r-o')
    grid on
    xlabel('SNR')
    ylabel('Average sum rate (bits per channel use)')
    set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1)
    title(sprintf('%s-IFC, K=%d, T=%d, R=%d, \\epsilon=1e-3', case_name, K, T, R),'Interpreter','tex');

    hold on
    plot(snr_range,rate_snr,'b-*')
    legend('WMMSE int10','WMMSE')
end
