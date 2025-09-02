function results = enhanced_dcc_sim(cfg)

%  Default configuration 
if nargin==0 || isempty(cfg)
    cfg.K_list    = [50 100 200 400];
    cfg.L         = 64;          % # APs
    cfg.N         = 4;           % antennas / AP
    cfg.tau_p     = 8;           % pilot length
    cfg.kappa_dB  = 5;           % Rician K‑factor [dB]
    cfg.cell      = 1000;        % area side [m]
    cfg.P_UL_dBm  = 15;
    cfg.P_DL_dBm  = 20;
    cfg.mobility  = "low";       % 'low' or 'high'
    cfg.gamma_min_dB = -3;
end

% Derived constants 
cfg.P_UL = 10^((cfg.P_UL_dBm-30)/10);   % W
cfg.P_DL = 10^((cfg.P_DL_dBm-30)/10);
Klist = cfg.K_list;  nK = numel(Klist);

%  Preallocate metrics 
SE_avg  = zeros(1,nK);
EE_avg  = zeros(1,nK);
Jain    = zeros(1,nK);
CPU_ms  = zeros(1,nK);
FH_mbps = zeros(1,nK);

%  MAIN LOOP over UE counts 
for ik = 1:nK
    K = Klist(ik);

    % (1) Generate random topology
    [ap_pos, ue_pos, beta_lin] = generate_topology(cfg, cfg.L, cfg.cell, K);

    % (2) Generate Rician channel matrices  (LOS+NLOS)
    Hhat = generate_rician_channel(cfg, ap_pos, ue_pos, beta_lin, K);

    % (3) Graph‑Coloring pilot + cluster formation
    [pilotIndex, Ck, Dl] = graph_coloring_pilot(cfg, beta_lin, K);

    % (4) MR‑stub UL/DL run  (replace with real LP‑MMSE + ELPC if available)
    [SE_k,EE_k,cpu_ms_iter,fh_bits_iter] = run_uplink_downlink(cfg,Hhat,Ck,Dl,K);

    % (5) Aggregate
    SE_avg(ik)  = mean(SE_k);
    EE_avg(ik)  = mean(EE_k);
    Jain(ik)    = (sum(SE_k))^2/(K*sum(SE_k.^2));
    CPU_ms(ik)  = cpu_ms_iter;      % already normalized per AP
    FH_mbps(ik) = fh_bits_iter/1e6;   % bits → Mbit/s

end

% Pack results 
results.cfg = cfg;
results.SE_avg = SE_avg;
results.EE_avg = EE_avg;
results.Jain = Jain;
results.cpu_ms = CPU_ms;
results.fh_mbps = FH_mbps;

%  Plot (fixed blank figure bug) 
fig = figure("Name","Enhanced DCC – Scalability & Performance");
layout = tiledlayout(2,2,"TileSpacing","Compact","Padding","Compact");

nexttile; plot(Klist,SE_avg,'-o','LineWidth',1.4); grid on;
xlabel('K'); ylabel('Avg SE (bit/s/Hz)'); title('Spectral Efficiency');

nexttile; plot(Klist,EE_avg,'-s','LineWidth',1.4); grid on;
xlabel('K'); ylabel('Avg EE (bit/J)'); title('Energy Efficiency');

nexttile; semilogx(Klist,CPU_ms,'-x','LineWidth',1.4); grid on;
xlabel('K'); ylabel('CPU (ms/AP)'); title('Computation');

nexttile; semilogx(Klist,FH_mbps,'-*','LineWidth',1.4); grid on;
xlabel('K'); ylabel('FH (Mbit/s/AP)'); title('Fronthaul');

title(layout, 'Enhanced DCC ‑ Scalability & Performance');
results.fig = fig;

end

%% generate random topology 
function [ap_pos, ue_pos, beta_lin] = generate_topology(cfg,L,cellSide,K)
% ap_pos[Lx2], ue_pos[Kx2] uniformly random; beta ~ d^{-3.76}
N0 = 1; pathExp = 3.76;
ap_pos = cellSide*rand(L,2);
ue_pos = cellSide*rand(K,2);
for l=1:L
   d = sqrt(sum((ue_pos - ap_pos(l,:)).^2,2))+1; % +1 to avoid 0
   beta_lin(:,l) = (N0./d).^pathExp;
end
end

%%  Rician channel (placeholder) 
function Hhat = generate_rician_channel(cfg,ap_pos,ue_pos,beta_lin,K)


L  = cfg.L;                 % #APs
N  = cfg.N;                 % antennas per AP
M  = L*N;                   % total antennas
kappa = 10^(cfg.kappa_dB/10);
LOS_ratio  = sqrt(kappa/(1+kappa));
NLOS_ratio = sqrt(1/(1+kappa));

% ULA steering vector function 
a = @(ang) exp(1j*(0:N-1).'*pi*sin(ang));  % N×1 phase steering

%  Allocate 
Hhat = zeros(M,K);

for l = 1:L
    antIdx = (l-1)*N + (1:N);   % antenna indices in stacked vector
    % 取出大尺度增益 β(k,l)
    gain = sqrt(beta_lin(:,l).');        % 1×K

    %  AoA 計算：AP→UE 連線的水平角 
    vec = ue_pos - ap_pos(l,:);          % K×2 vector difference
    AoA = atan2(vec(:,2), vec(:,1));     % K×1  ∈ [-π,π]

    %  LOS component 
    Vlos = zeros(N,K);
    for k = 1:K
        Vlos(:,k) = a(AoA(k));           % N×1 steering
    end

    % NLOS component
    Gnlos = (randn(N,K) + 1j*randn(N,K))/sqrt(2);

    %  Combine & scale 
    Hhat(antIdx,:) = ( LOS_ratio  * Vlos + NLOS_ratio * Gnlos ) .* gain;
end
end
%% Helper: GCPA pilot (greedy DSATUR)  
function [pilotIndex,Ck,Dl] = graph_coloring_pilot(cfg,beta_lin,K)
tau = cfg.tau_p;   L = cfg.L;

%   UE-UE interference graph (βil·βjl > ξ)  
xi = 1e-11;                                  % 門檻可自行調
G  = (beta_lin*beta_lin.') > xi;             % K×K boolean
G(1:K+1:end) = false;                        % 去自環

%  DSATUR 著色  
pilotIndex = zeros(1,K);
satDeg = zeros(1,K);                         % saturation degree
deg     = sum(G,2).';                        % 原始度數
for t = 1:K
    % 1) 找尚未上色中 saturation-deg 最大、再比度數
    uncol = find(pilotIndex==0);
    metric   = satDeg(uncol)*1e4 + deg(uncol);   % *1e4 保留排序優先權
    [~,pos]  = max(metric);                      % 只有一個輸入陣列
    k        = uncol(pos);                       % 真正的 UE index


    % 2) 給它最小可用顏色
    neighColor = pilotIndex(G(k,:) & pilotIndex>0);
    c = find(~ismember(1:tau,neighColor),1);
    if isempty(c), c = randi(tau); end       % 若超出，隨機給
    pilotIndex(k)=c;
    % 3) 更新 saturation 度
    neigh = find(G(k,:));
    for n = neigh
        if pilotIndex(n)==0
            satDeg(n)=numel(unique([pilotIndex(G(n,:) & pilotIndex>0),c]));
        end
    end
end

%  user-centric cluster: 每 UE 取前三大 β 
[~,idxSort] = sort(beta_lin,2,'descend');
Ck = arrayfun(@(k) idxSort(k,1:min(3,L)), 1:K,'uni',0);

% AP-centric list Dl 
Dl = arrayfun(@(l) find(cellfun(@(c) any(c==l), Ck)), 1:L,'uni',0);

% Enforce |Dl(l)| ≤ τp 
tau = cfg.tau_p;
for l = 1:L
    while numel(Dl{l}) > tau
        users = Dl{l};                       % UE list served by AP l
        % 踢掉 β 最小的那個 UE (最不影響效能)
        [~,idxMin] = min(beta_lin(users,l));
        uKick = users(idxMin);
        % 從 UE 的 cluster 移除 AP l
        Ck{uKick}(Ck{uKick}==l) = [];
        % 從 Dl{l} 移除該 UE
        Dl{l}(idxMin) = [];
    end
end



end


%% UL/DL

function [SE_k,EE_k,cpu_ms,fh_bits] = run_uplink_downlink(cfg,Hhat,Ck,Dl,K)
L = cfg.L; N = cfg.N; P = cfg.P_UL;
SE_k = zeros(1,K); cpu_each = zeros(1,L);

for l = 1:L
    tic
    users = Dl{l};             % ≤ τp
    if isempty(users), cpu_each(l)=0; continue; end
    ant = (l-1)*N + (1:N);
    H  = Hhat(ant,users);      % N × |users|
    V  = H;                    % MR
    sig = P*abs(diag(V'*H)).'.^2;
    int = P*(sum(abs(V'*H).^2,2).' - sig);
    noise = sum(abs(V).^2,1)*1e-13;
    SINR = sig./(int+noise);
    SE_k(users) = (1-cfg.tau_p/200).*log2(1+SINR);
    cpu_each(l) = toc*1000;    % ms
end

EE_k  = SE_k./(P+0.1);
cpu_ms = mean(cpu_each);       % ≈ 常數
fh_bits = 2*cfg.tau_p*32*1000; % 常數 512 kbit/s
end

