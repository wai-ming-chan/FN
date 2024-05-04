clc; clear all; close all;
% load fh.mat
rng(12348)

N_realization = 800;
SNRdB_list = [-2 0 2 4 6, 6.2:0.2:7];
% SNRdB_list = [2 6 10 16 20];
% SNRdB_list = [7];


Ns = 1620; % total number of sc
Nh = 14;    % number of OFDM symbol in one slot time
Nu1 = 32; % number of URLLC UE
Nu2 = 3000; % number of mMTC UE
Nu = Nu1 + Nu2;
Nj = 32*16; % number of jammed sc
Nj_band = 512; % number of sc in jamming band
Ni1 = 1; % demanded number of sc for URLLC UE
Ni2 = 1; % demanded number of sc for mMTC UE

Xi_noFH = genX_noFH(Ns,Nh,Nu);

% J = zeros(Ns, Nh);
% % =============================== Uniform Jamming ===============================
% mid_jb = ceil(Ns/2); 
% Nj_lb = mid_jb - Nj_band /2;
% Nj_ub = mid_jb + Nj_band /2;

% for n=1:Nh
%   jammed_sc_idx = Nj_lb + randperm(Nj_band, Nj);
%   J(jammed_sc_idx,n) = 1;
% end

% =============================== Gauss Jamming ===============================
pj = zeros(Ns,1);
pi_mpf_count = 0;

mid_jb = ceil(Ns/2); 
mu_jam = mid_jb;
sig_jam = Nj_band/ (4); % [1/2, 1/4, 1/6] -> [0.68, 0.95, 99.7] interval

pj_obj.a = 1;
pj_obj.b = Ns;
pj_obj.mu = mu_jam;
pj_obj.sigma = sig_jam;
pj_obj.n = Nj;

pjs=zeros(Ns,N_realization);
pj = zeros(Ns,1);
for i_real = 1:N_realization+2
    Jset=[];
    
    for n=1:Nh
        jam_pattern = genJamPattern(pj_obj);
        Jset = [Jset, jam_pattern];
    end
    Jsets{i_real} = Jset;

    pj_tmp = genPMF_pj(Jset(:), pj_obj.b);
    pj = (i_real-1)/i_real * pj + (1/i_real) * pj_tmp;
    
    pjs(:,i_real) = pj;
end
clear pj_tmp pj jam_pattern Jset 
pj_true = genTrue_pj(pj_obj)';

% ===============================  ===============================

for i_real = 1:N_realization
    fprintf('realization: %d \t.\n', i_real);
    Jset = Jsets{i_real+1};
    pj = pjs(:,i_real);
    %pj = pj_true;

%     N_jam_training = 14;
%     jam_samples=[];
%     for i=1:N_jam_training
%         jam_pattern = genJamPattern(pj_obj);
%         jam_samples = [jam_samples; jam_pattern];
%     end
%     pj = genPMF_pj( jam_samples, pj_obj.b );

  % ================ (1,2) Full Knowledge of Jamming FH pattern case ===============
  Fset = 1:Ns;

  Xi =zeros(Nu,Nh);
  UEoffset = 1;
  for n=1:Nh
    JIset_n = Jset(:,n)';
    for ii=1:Nu
      i = mod(ii-1 + UEoffset-1, Nu)+1;
      FJIset = setdiff(Fset, JIset_n);

      if (length(FJIset) > 0) % not equal to f_null
        Xi(i,n) = FJIset(randi(length(FJIset), 1,1) );
        JIset_n = [JIset_n, Xi(i,n)];
        UEoffset_tmp = i;
      else
        Xi(i,n) = 0; % null
      end
      if ii == Nu
        UEoffset = UEoffset_tmp;
      end
    end
  end
  Xi_Det = Xi;
  
  % fprintf('End Deterministic case\n');

  % ================ Statisctical Jamming FH pattern case =============== 
  % ================ (3) Without Prior of jamming distribution =======================
  Xi =zeros(Nu,Nh);
  
  
  UEoffset = 1;
  for n=1:Nh
    p_JI = pj; % initial condition
    
    ue_permuted = randperm(Ns-Nj);
    i_prev = 0;

    yi_prime = zeros(Ns,1);
    for ii=1:Ns-Nj
      i = mod( ue_permuted(ii)-1 + UEoffset-1, Nu) + 1;
      
      oneHot_i1 = zeros(Ns,1);
      if (ii >=2 && Xi(i_prev,n) ~= 0)
        oneHot_i1(Xi(i_prev,n), 1) = 1;
      end 
      % p_JgivenI = pj + (1 - pj).* oneHot_i1;
      
      if ii > 2
         % p_JI = (ii-2)/(ii-1) * p_JI + 1/(ii-1) * p_JgivenI;
         p_JI = (Nj + ii-2)/(Nj+ii-1) * p_JI + 1/(Nj+ii-1) * oneHot_i1;
      elseif ii == 2
        p_JI = (Nj / Nj+1) * pj + 1/(Nj+1) * oneHot_i1;
      else
        p_JI = pj;
      end

      [sortVal, sortIdx] = sort(p_JI); 
      avail_idx = sortIdx( find(sortVal == sortVal(1)) );
      uei_sc_idx = randperm(length(avail_idx), 1);
      Xi(i,n) = (avail_idx(uei_sc_idx));
     
      % [~, minIdx] = min(p_JI);
      % Xi(i,n) = (minIdx);
      i_prev = i;
    end
    
    UEoffset = mod(UEoffset+ Ns-Nj-1, Nu)+1;
  end
  Xi_Stat = Xi;



  % ================ (4) With Prior of jamming distribution =======================
  Xi_m2 =zeros(Nu,Nh);

  UEoffset = 1;
  for n=1:Nh
    p_JI = pj_true; % initial condition
    
    ue_permuted = randperm(Ns-Nj);
    i_prev = 0;

    yi_prime = zeros(Ns,1);
    for ii=1:Ns-Nj
      i = mod( ue_permuted(ii)-1 + UEoffset-1, Nu) + 1;
      
      oneHot_i1 = zeros(Ns,1);
      if (ii >=2 && Xi_m2(i_prev,n) ~= 0)
        oneHot_i1(Xi_m2(i_prev,n), 1) = 1;
      end 
      % p_JgivenI = pj + (1 - pj).* oneHot_i1;
      
      if ii > 2
         % p_JI = (ii-2)/(ii-1) * p_JI + 1/(ii-1) * p_JgivenI;
         p_JI = (Nj + ii-2)/(Nj+ii-1) * p_JI + 1/(Nj+ii-1) * oneHot_i1;
      elseif ii == 2
        p_JI = (Nj / Nj+1) * pj_true + 1/(Nj+1) * oneHot_i1;
      else
        p_JI = pj_true;
      end

      [sortVal, sortIdx] = sort(p_JI); 
      avail_idx = sortIdx( find(sortVal == sortVal(1)) );
      uei_sc_idx = randperm(length(avail_idx), 1);
      Xi_m2(i,n) = (avail_idx(uei_sc_idx));
     
      % [~, minIdx] = min(p_JI);
      % Xi(i,n) = (minIdx);
      i_prev = i;
    end
    
    UEoffset = mod(UEoffset+ Ns-Nj-1, Nu)+1;
  end
  
  Xi_Stat_prior = Xi_m2; 
 
  % fprintf('End statistic case\n');
  % ================ (5) Uniform FH Patterns =======================
  Xi =zeros(Nu,Nh);
  
  Fset = 1:Ns;

  Xi =zeros(Nu,Nh);
  UEoffset = 1;
  for n=1:Nh
    JIset_n = [];
    for ii=1:Nu
      i = mod(ii-1 + UEoffset-1, Nu)+1;
      FJIset = setdiff(Fset, JIset_n);

      if (length(FJIset) > 0) % not equal to f_null
        Xi(i,n) = FJIset(randi(length(FJIset), 1,1) );
        JIset_n = [JIset_n, Xi(i,n)];
        UEoffset_tmp = i;
      else
        Xi(i,n) = 0; % null
      end
      if ii == Nu
        UEoffset = UEoffset_tmp;
      end
    end
  end
  
  Xi_Uniform = Xi; 
 
  % fprintf('End statistic case\n');


  % ================ Non-causal case =====================================
  % X_jam = Jset(:,1:Nh);
  X_jam = Jsets{i_real+1};
  % oo = transX2Y(Xi_Stat, X_jam, Ns);
  oo = transX2Y(Xi_Det, X_jam, Ns);
  Yi = oo.Yi;
  Y_jam = oo.Y_jam;

  %SNRdB_list = [-2 0 2 4 6 7];
  
  for i_snr = 1:length(SNRdB_list)
    SNRdB = SNRdB_list(i_snr);
    output = simQPSK_Y(Yi, Y_jam, SNRdB);
    nz_SINR = output.SINR_ues;
    
    
    y1_SINR_min(i_snr) = min(nz_SINR);        y1_SINR_mean(i_snr) = mean(nz_SINR);
    y1_SINR_p995(i_snr) = prctile(nz_SINR,99.5); 
    y1_SINR_p99(i_snr) = prctile(nz_SINR,99); y1_SINR_p95(i_snr) = prctile(nz_SINR,95);
    y1_SINR_p90(i_snr) = prctile(nz_SINR,90); y1_SINR_p80(i_snr) = prctile(nz_SINR,80);
    y1_SINR_p70(i_snr) = prctile(nz_SINR,70); y1_SINR_p60(i_snr) = prctile(nz_SINR,60);
    y1_SINR_p50(i_snr) = prctile(nz_SINR,50); 
    y1_SINR_p40(i_snr) = prctile(nz_SINR,40); y1_SINR_p30(i_snr) = prctile(nz_SINR,30);
    y1_SINR_p20(i_snr) = prctile(nz_SINR,20); y1_SINR_p10(i_snr) = prctile(nz_SINR,10);
    y1_SINR_p05(i_snr) = prctile(nz_SINR, 5); y1_SINR_p01(i_snr) = prctile(nz_SINR, 1); 
    y1_SINR_p005(i_snr) = prctile(nz_SINR,0.5);
    
  
    output = simBER_givenYiYjam(Yi, Y_jam, SNRdB);
    nz_BER = output.BER;
    y1_BER_min(i_snr) = min(nz_BER);        y1_BER_max(i_snr) = max(nz_BER);
    y1_BER_p995(i_snr) = prctile(nz_BER,99.5);
    y1_BER_p99(i_snr) = prctile(nz_BER,99); y1_BER_p95(i_snr) = prctile(nz_BER,95);
    y1_BER_p90(i_snr) = prctile(nz_BER,90); y1_BER_p80(i_snr) = prctile(nz_BER,80);
    y1_BER_p70(i_snr) = prctile(nz_BER,70); y1_BER_p60(i_snr) = prctile(nz_BER,60);
    y1_BER_p50(i_snr) = prctile(nz_BER,50); 
    y1_BER_p40(i_snr) = prctile(nz_BER,40); y1_BER_p30(i_snr) = prctile(nz_BER,30);
    y1_BER_p20(i_snr) = prctile(nz_BER,20); y1_BER_p10(i_snr) = prctile(nz_BER,10);
    y1_BER_p05(i_snr) = prctile(nz_BER, 5); y1_BER_p01(i_snr) = prctile(nz_BER, 1);
    y1_BER_p005(i_snr) = prctile(nz_BER,0.5);
    

  end
  yDet_noncausal_SINR_min(i_real, :) = mean(y1_SINR_min,1); 
  yDet_noncausal_SINR_mean(i_real,:) = mean(y1_SINR_mean,1);
  yDet_noncausal_SINR_p995(i_real, :) = mean(y1_SINR_p995,1); 
  yDet_noncausal_SINR_p99(i_real, :) = mean(y1_SINR_p99,1); yDet_noncausal_SINR_p95(i_real, :) = mean(y1_SINR_p95,1); 
  yDet_noncausal_SINR_p90(i_real, :) = mean(y1_SINR_p90,1); yDet_noncausal_SINR_p80(i_real, :) = mean(y1_SINR_p80,1);
  yDet_noncausal_SINR_p70(i_real, :) = mean(y1_SINR_p70,1); yDet_noncausal_SINR_p60(i_real, :) = mean(y1_SINR_p60,1); 
  yDet_noncausal_SINR_p50(i_real, :) = mean(y1_SINR_p50,1);
  yDet_noncausal_SINR_p40(i_real, :) = mean(y1_SINR_p40,1); yDet_noncausal_SINR_p30(i_real, :) = mean(y1_SINR_p30,1);
  yDet_noncausal_SINR_p20(i_real, :) = mean(y1_SINR_p20,1); yDet_noncausal_SINR_p10(i_real, :) = mean(y1_SINR_p10,1);
  yDet_noncausal_SINR_p05(i_real, :) = mean(y1_SINR_p05,1); yDet_noncausal_SINR_p01(i_real, :) = mean(y1_SINR_p01,1);
  yDet_noncausal_SINR_p005(i_real, :) = mean(y1_SINR_p005,1);
  
  
  yDet_noncausal_BER_min(i_real, :) = mean(y1_BER_min,1);
  yDet_noncausal_BER_max(i_real, :) = mean(y1_BER_max,1);
  yDet_noncausal_BER_p995(i_real, :) = mean(y1_BER_p995,1);
  yDet_noncausal_BER_p99(i_real, :) = mean(y1_BER_p99,1); yDet_noncausal_BER_p95(i_real, :) = mean(y1_BER_p95,1);
  yDet_noncausal_BER_p90(i_real, :) = mean(y1_BER_p90,1); yDet_noncausal_BER_p80(i_real, :) = mean(y1_BER_p80,1);
  yDet_noncausal_BER_p70(i_real, :) = mean(y1_BER_p70,1); yDet_noncausal_BER_p60(i_real, :) = mean(y1_BER_p60,1);
  yDet_noncausal_BER_p50(i_real, :) = mean(y1_BER_p50,1);
  yDet_noncausal_BER_p40(i_real, :) = mean(y1_BER_p40,1); yDet_noncausal_BER_p30(i_real, :) = mean(y1_BER_p30,1);
  yDet_noncausal_BER_p20(i_real, :) = mean(y1_BER_p20,1); yDet_noncausal_BER_p10(i_real, :) = mean(y1_BER_p10,1);
  yDet_noncausal_BER_p05(i_real, :) = mean(y1_BER_p05,1); yDet_noncausal_BER_p01(i_real, :) = mean(y1_BER_p01,1);
  yDet_noncausal_BER_p005(i_real, :) = mean(y1_BER_p005,1);
  
  % ================ Causal case =========================================
%   X_jam = Jset(:,2:Nh+1);
  X_jam = Jsets{i_real+2};
  oo = transX2Y(Xi_Det, X_jam, Ns);
  Yi = oo.Yi;
  Y_jam = oo.Y_jam;

%   SNRdB_list = [-2 0 2 4 6 7];
  for i_snr = 1:length(SNRdB_list)
    SNRdB = SNRdB_list(i_snr);
    output = simQPSK_Y(Yi, Y_jam, SNRdB);
    nz_SINR = output.SINR_ues;

    y1_SINR_min(i_snr) = min(nz_SINR);        y1_SINR_mean(i_snr) = mean(nz_SINR);
    y1_SINR_p995(i_snr) = prctile(nz_SINR,99.5); 
    y1_SINR_p99(i_snr) = prctile(nz_SINR,99); y1_SINR_p95(i_snr) = prctile(nz_SINR,95);
    y1_SINR_p90(i_snr) = prctile(nz_SINR,90); y1_SINR_p80(i_snr) = prctile(nz_SINR,80);
    y1_SINR_p70(i_snr) = prctile(nz_SINR,70); y1_SINR_p60(i_snr) = prctile(nz_SINR,60);
    y1_SINR_p50(i_snr) = prctile(nz_SINR,50); 
    y1_SINR_p40(i_snr) = prctile(nz_SINR,40); y1_SINR_p30(i_snr) = prctile(nz_SINR,30);
    y1_SINR_p20(i_snr) = prctile(nz_SINR,20); y1_SINR_p10(i_snr) = prctile(nz_SINR,10);
    y1_SINR_p05(i_snr) = prctile(nz_SINR, 5); y1_SINR_p01(i_snr) = prctile(nz_SINR, 1); 
    y1_SINR_p005(i_snr) = prctile(nz_SINR,0.5);


    output = simBER_givenYiYjam(Yi, Y_jam, SNRdB);
    nz_BER = output.BER;
    y1_BER_min(i_snr) = min(nz_BER);        y1_BER_max(i_snr) = max(nz_BER);
    y1_BER_p995(i_snr) = prctile(nz_BER,99.5);
    y1_BER_p99(i_snr) = prctile(nz_BER,99); y1_BER_p95(i_snr) = prctile(nz_BER,95);
    y1_BER_p90(i_snr) = prctile(nz_BER,90); y1_BER_p80(i_snr) = prctile(nz_BER,80);
    y1_BER_p70(i_snr) = prctile(nz_BER,70); y1_BER_p60(i_snr) = prctile(nz_BER,60);
    y1_BER_p50(i_snr) = prctile(nz_BER,50); 
    y1_BER_p40(i_snr) = prctile(nz_BER,40); y1_BER_p30(i_snr) = prctile(nz_BER,30);
    y1_BER_p20(i_snr) = prctile(nz_BER,20); y1_BER_p10(i_snr) = prctile(nz_BER,10);
    y1_BER_p05(i_snr) = prctile(nz_BER, 5); y1_BER_p01(i_snr) = prctile(nz_BER, 1);
    y1_BER_p005(i_snr) = prctile(nz_BER,0.5);

  end 
  yDet_causal_SINR_min(i_real, :) = mean(y1_SINR_min,1); 
  yDet_causal_SINR_mean(i_real,:) = mean(y1_SINR_mean,1);
  yDet_causal_SINR_p995(i_real, :) = mean(y1_SINR_p995,1); 
  yDet_causal_SINR_p99(i_real, :) = mean(y1_SINR_p99,1); yDet_causal_SINR_p95(i_real, :) = mean(y1_SINR_p95,1); 
  yDet_causal_SINR_p90(i_real, :) = mean(y1_SINR_p90,1); yDet_causal_SINR_p80(i_real, :) = mean(y1_SINR_p80,1);
  yDet_causal_SINR_p70(i_real, :) = mean(y1_SINR_p70,1); yDet_causal_SINR_p60(i_real, :) = mean(y1_SINR_p60,1); 
  yDet_causal_SINR_p50(i_real, :) = mean(y1_SINR_p50,1);
  yDet_causal_SINR_p40(i_real, :) = mean(y1_SINR_p40,1); yDet_causal_SINR_p30(i_real, :) = mean(y1_SINR_p30,1);
  yDet_causal_SINR_p20(i_real, :) = mean(y1_SINR_p20,1); yDet_causal_SINR_p10(i_real, :) = mean(y1_SINR_p10,1);
  yDet_causal_SINR_p05(i_real, :) = mean(y1_SINR_p05,1); yDet_causal_SINR_p01(i_real, :) = mean(y1_SINR_p01,1);
  yDet_causal_SINR_p005(i_real, :) = mean(y1_SINR_p005,1);

  
  
  yDet_causal_BER_min(i_real, :) = mean(y1_BER_min,1);
  yDet_causal_BER_max(i_real, :) = mean(y1_BER_max,1);
  yDet_causal_BER_p995(i_real, :) = mean(y1_BER_p995,1);
  yDet_causal_BER_p99(i_real, :) = mean(y1_BER_p99,1); yDet_causal_BER_p95(i_real, :) = mean(y1_BER_p95,1);
  yDet_causal_BER_p90(i_real, :) = mean(y1_BER_p90,1); yDet_causal_BER_p80(i_real, :) = mean(y1_BER_p80,1);
  yDet_causal_BER_p70(i_real, :) = mean(y1_BER_p70,1); yDet_causal_BER_p60(i_real, :) = mean(y1_BER_p60,1);
  yDet_causal_BER_p50(i_real, :) = mean(y1_BER_p50,1);
  yDet_causal_BER_p40(i_real, :) = mean(y1_BER_p40,1); yDet_causal_BER_p30(i_real, :) = mean(y1_BER_p30,1);
  yDet_causal_BER_p20(i_real, :) = mean(y1_BER_p20,1); yDet_causal_BER_p10(i_real, :) = mean(y1_BER_p10,1);
  yDet_causal_BER_p05(i_real, :) = mean(y1_BER_p05,1); yDet_causal_BER_p01(i_real, :) = mean(y1_BER_p01,1);
  yDet_causal_BER_p005(i_real, :) = mean(y1_BER_p005,1);

  % ================ Statistical case =====================================
  % X_jam = Jset(:,1:Nh);
  X_jam = Jsets{i_real+1};
  oo = transX2Y(Xi_Stat, X_jam, Ns);
  Yi = oo.Yi;
  Y_jam = oo.Y_jam;

%   SNRdB_list = [-2 0 2 4 6 7];
  for i_snr = 1:length(SNRdB_list)
    SNRdB = SNRdB_list(i_snr);
    output = simQPSK_Y(Yi, Y_jam, SNRdB);
    nz_SINR = output.SINR_ues;
    
    y1_SINR_min(i_snr) = min(nz_SINR);        y1_SINR_mean(i_snr) = mean(nz_SINR);
    y1_SINR_p995(i_snr) = prctile(nz_SINR,99.5); 
    y1_SINR_p99(i_snr) = prctile(nz_SINR,99); y1_SINR_p95(i_snr) = prctile(nz_SINR,95);
    y1_SINR_p90(i_snr) = prctile(nz_SINR,90); y1_SINR_p80(i_snr) = prctile(nz_SINR,80);
    y1_SINR_p70(i_snr) = prctile(nz_SINR,70); y1_SINR_p60(i_snr) = prctile(nz_SINR,60);
    y1_SINR_p50(i_snr) = prctile(nz_SINR,50); 
    y1_SINR_p40(i_snr) = prctile(nz_SINR,40); y1_SINR_p30(i_snr) = prctile(nz_SINR,30);
    y1_SINR_p20(i_snr) = prctile(nz_SINR,20); y1_SINR_p10(i_snr) = prctile(nz_SINR,10);
    y1_SINR_p05(i_snr) = prctile(nz_SINR, 5); y1_SINR_p01(i_snr) = prctile(nz_SINR, 1); 
    y1_SINR_p005(i_snr) = prctile(nz_SINR,0.5);

    output = simBER_givenYiYjam(Yi, Y_jam, SNRdB);
    nz_BER = output.BER;
    y1_BER_min(i_snr) = min(nz_BER);        y1_BER_max(i_snr) = max(nz_BER);
    y1_BER_p995(i_snr) = prctile(nz_BER,99.5);
    y1_BER_p99(i_snr) = prctile(nz_BER,99); y1_BER_p95(i_snr) = prctile(nz_BER,95);
    y1_BER_p90(i_snr) = prctile(nz_BER,90); y1_BER_p80(i_snr) = prctile(nz_BER,80);
    y1_BER_p70(i_snr) = prctile(nz_BER,70); y1_BER_p60(i_snr) = prctile(nz_BER,60);
    y1_BER_p50(i_snr) = prctile(nz_BER,50); 
    y1_BER_p40(i_snr) = prctile(nz_BER,40); y1_BER_p30(i_snr) = prctile(nz_BER,30);
    y1_BER_p20(i_snr) = prctile(nz_BER,20); y1_BER_p10(i_snr) = prctile(nz_BER,10);
    y1_BER_p05(i_snr) = prctile(nz_BER, 5); y1_BER_p01(i_snr) = prctile(nz_BER, 1);
    y1_BER_p005(i_snr) = prctile(nz_BER,0.5);
  end 

  yStat_SINR_min(i_real, :) = mean(y1_SINR_min,1); 
  yStat_SINR_mean(i_real,:) = mean(y1_SINR_mean,1);
  yStat_SINR_p995(i_real, :) = mean(y1_SINR_p995,1); 
  yStat_SINR_p99(i_real, :) = mean(y1_SINR_p99,1); yStat_SINR_p95(i_real, :) = mean(y1_SINR_p95,1); 
  yStat_SINR_p90(i_real, :) = mean(y1_SINR_p90,1); yStat_SINR_p80(i_real, :) = mean(y1_SINR_p80,1);
  yStat_SINR_p70(i_real, :) = mean(y1_SINR_p70,1); yStat_SINR_p60(i_real, :) = mean(y1_SINR_p60,1); 
  yStat_SINR_p50(i_real, :) = mean(y1_SINR_p50,1);
  yStat_SINR_p40(i_real, :) = mean(y1_SINR_p40,1); yStat_SINR_p30(i_real, :) = mean(y1_SINR_p30,1);
  yStat_SINR_p20(i_real, :) = mean(y1_SINR_p20,1); yStat_SINR_p10(i_real, :) = mean(y1_SINR_p10,1);
  yStat_SINR_p05(i_real, :) = mean(y1_SINR_p05,1); yStat_SINR_p01(i_real, :) = mean(y1_SINR_p01,1);
  yStat_SINR_p005(i_real, :) = mean(y1_SINR_p005,1);

  yStat_BER_min(i_real, :) = mean(y1_BER_min,1);
  yStat_BER_max(i_real, :) = mean(y1_BER_max,1);
  yStat_BER_p995(i_real, :) = mean(y1_BER_p995,1);
  yStat_BER_p99(i_real, :) = mean(y1_BER_p99,1); yStat_BER_p95(i_real, :) = mean(y1_BER_p95,1);
  yStat_BER_p90(i_real, :) = mean(y1_BER_p90,1); yStat_BER_p80(i_real, :) = mean(y1_BER_p80,1);
  yStat_BER_p70(i_real, :) = mean(y1_BER_p70,1); yStat_BER_p60(i_real, :) = mean(y1_BER_p60,1);
  yStat_BER_p50(i_real, :) = mean(y1_BER_p50,1);
  yStat_BER_p40(i_real, :) = mean(y1_BER_p40,1); yStat_BER_p30(i_real, :) = mean(y1_BER_p30,1);
  yStat_BER_p20(i_real, :) = mean(y1_BER_p20,1); yStat_BER_p10(i_real, :) = mean(y1_BER_p10,1);
  yStat_BER_p05(i_real, :) = mean(y1_BER_p05,1); yStat_BER_p01(i_real, :) = mean(y1_BER_p01,1);
  yStat_BER_p005(i_real, :) = mean(y1_BER_p005,1);

  

  % ================ Statistical case (with prior) =====================================
  %X_jam = Jset(:,1:Nh);
  X_jam = Jsets{i_real+1};
  oo = transX2Y(Xi_Stat_prior, X_jam, Ns);
  Yi = oo.Yi;
  Y_jam = oo.Y_jam;

%   SNRdB_list = [-2 0 2 4 6 7];
  for i_snr = 1:length(SNRdB_list)
    SNRdB = SNRdB_list(i_snr);
    output = simQPSK_Y(Yi, Y_jam, SNRdB);
    nz_SINR = output.SINR_ues;
    
    y1_SINR_min(i_snr) = min(nz_SINR);        y1_SINR_mean(i_snr) = mean(nz_SINR);
    y1_SINR_p995(i_snr) = prctile(nz_SINR,99.5); 
    y1_SINR_p99(i_snr) = prctile(nz_SINR,99); y1_SINR_p95(i_snr) = prctile(nz_SINR,95);
    y1_SINR_p90(i_snr) = prctile(nz_SINR,90); y1_SINR_p80(i_snr) = prctile(nz_SINR,80);
    y1_SINR_p70(i_snr) = prctile(nz_SINR,70); y1_SINR_p60(i_snr) = prctile(nz_SINR,60);
    y1_SINR_p50(i_snr) = prctile(nz_SINR,50); 
    y1_SINR_p40(i_snr) = prctile(nz_SINR,40); y1_SINR_p30(i_snr) = prctile(nz_SINR,30);
    y1_SINR_p20(i_snr) = prctile(nz_SINR,20); y1_SINR_p10(i_snr) = prctile(nz_SINR,10);
    y1_SINR_p05(i_snr) = prctile(nz_SINR, 5); y1_SINR_p01(i_snr) = prctile(nz_SINR, 1); 
    y1_SINR_p005(i_snr) = prctile(nz_SINR,0.5);

    output = simBER_givenYiYjam(Yi, Y_jam, SNRdB);
    nz_BER = output.BER;
    y1_BER_min(i_snr) = min(nz_BER);        y1_BER_max(i_snr) = max(nz_BER);
    y1_BER_p995(i_snr) = prctile(nz_BER,99.5);
    y1_BER_p99(i_snr) = prctile(nz_BER,99); y1_BER_p95(i_snr) = prctile(nz_BER,95);
    y1_BER_p90(i_snr) = prctile(nz_BER,90); y1_BER_p80(i_snr) = prctile(nz_BER,80);
    y1_BER_p70(i_snr) = prctile(nz_BER,70); y1_BER_p60(i_snr) = prctile(nz_BER,60);
    y1_BER_p50(i_snr) = prctile(nz_BER,50); 
    y1_BER_p40(i_snr) = prctile(nz_BER,40); y1_BER_p30(i_snr) = prctile(nz_BER,30);
    y1_BER_p20(i_snr) = prctile(nz_BER,20); y1_BER_p10(i_snr) = prctile(nz_BER,10);
    y1_BER_p05(i_snr) = prctile(nz_BER, 5); y1_BER_p01(i_snr) = prctile(nz_BER, 1);
    y1_BER_p005(i_snr) = prctile(nz_BER,0.5);
    

  end 
  yStat_m2_SINR_min(i_real, :) = mean(y1_SINR_min,1); 
  yStat_m2_SINR_mean(i_real,:) = mean(y1_SINR_mean,1);
  yStat_m2_SINR_p995(i_real, :) = mean(y1_SINR_p995,1); 
  yStat_m2_SINR_p99(i_real, :) = mean(y1_SINR_p99,1); yStat_m2_SINR_p95(i_real, :) = mean(y1_SINR_p95,1); 
  yStat_m2_SINR_p90(i_real, :) = mean(y1_SINR_p90,1); yStat_m2_SINR_p80(i_real, :) = mean(y1_SINR_p80,1);
  yStat_m2_SINR_p70(i_real, :) = mean(y1_SINR_p70,1); yStat_m2_SINR_p60(i_real, :) = mean(y1_SINR_p60,1); 
  yStat_m2_SINR_p50(i_real, :) = mean(y1_SINR_p50,1);
  yStat_m2_SINR_p40(i_real, :) = mean(y1_SINR_p40,1); yStat_m2_SINR_p30(i_real, :) = mean(y1_SINR_p30,1);
  yStat_m2_SINR_p20(i_real, :) = mean(y1_SINR_p20,1); yStat_m2_SINR_p10(i_real, :) = mean(y1_SINR_p10,1);
  yStat_m2_SINR_p05(i_real, :) = mean(y1_SINR_p05,1); yStat_m2_SINR_p01(i_real, :) = mean(y1_SINR_p01,1);
  yStat_m2_SINR_p005(i_real, :) = mean(y1_SINR_p005,1);

  yStat_m2_BER_min(i_real, :) = mean(y1_BER_min,1);
  yStat_m2_BER_max(i_real, :) = mean(y1_BER_max,1);
  yStat_m2_BER_p995(i_real, :) = mean(y1_BER_p995,1);
  yStat_m2_BER_p99(i_real, :) = mean(y1_BER_p99,1); yStat_m2_BER_p95(i_real, :) = mean(y1_BER_p95,1);
  yStat_m2_BER_p90(i_real, :) = mean(y1_BER_p90,1); yStat_m2_BER_p80(i_real, :) = mean(y1_BER_p80,1);
  yStat_m2_BER_p70(i_real, :) = mean(y1_BER_p70,1); yStat_m2_BER_p60(i_real, :) = mean(y1_BER_p60,1);
  yStat_m2_BER_p50(i_real, :) = mean(y1_BER_p50,1);
  yStat_m2_BER_p40(i_real, :) = mean(y1_BER_p40,1); yStat_m2_BER_p30(i_real, :) = mean(y1_BER_p30,1);
  yStat_m2_BER_p20(i_real, :) = mean(y1_BER_p20,1); yStat_m2_BER_p10(i_real, :) = mean(y1_BER_p10,1);
  yStat_m2_BER_p05(i_real, :) = mean(y1_BER_p05,1); yStat_m2_BER_p01(i_real, :) = mean(y1_BER_p01,1);
  yStat_m2_BER_p005(i_real, :) = mean(y1_BER_p005,1);

  
  
  % ================ Uniform FH case =====================================
  % X_jam = Jset(:,1:Nh);
  X_jam = Jsets{i_real+1};
  % oo = transX2Y(Xi_Stat, X_jam, Ns);
  oo = transX2Y(Xi_Uniform, X_jam, Ns);
  Yi = oo.Yi;
  Y_jam = oo.Y_jam;

%   SNRdB_list = [-2 0 2 4 6 7];
  for i_snr = 1:length(SNRdB_list)
    SNRdB = SNRdB_list(i_snr);
    output = simQPSK_Y(Yi, Y_jam, SNRdB);
    nz_SINR = output.SINR_ues;
    
    
    y1_SINR_min(i_snr) = min(nz_SINR);        y1_SINR_mean(i_snr) = mean(nz_SINR);
    y1_SINR_p995(i_snr) = prctile(nz_SINR,99.5); 
    y1_SINR_p99(i_snr) = prctile(nz_SINR,99); y1_SINR_p95(i_snr) = prctile(nz_SINR,95);
    y1_SINR_p90(i_snr) = prctile(nz_SINR,90); y1_SINR_p80(i_snr) = prctile(nz_SINR,80);
    y1_SINR_p70(i_snr) = prctile(nz_SINR,70); y1_SINR_p60(i_snr) = prctile(nz_SINR,60);
    y1_SINR_p50(i_snr) = prctile(nz_SINR,50); 
    y1_SINR_p40(i_snr) = prctile(nz_SINR,40); y1_SINR_p30(i_snr) = prctile(nz_SINR,30);
    y1_SINR_p20(i_snr) = prctile(nz_SINR,20); y1_SINR_p10(i_snr) = prctile(nz_SINR,10);
    y1_SINR_p05(i_snr) = prctile(nz_SINR, 5); y1_SINR_p01(i_snr) = prctile(nz_SINR, 1); 
    y1_SINR_p005(i_snr) = prctile(nz_SINR,0.5);

    output = simBER_givenYiYjam(Yi, Y_jam, SNRdB);
    nz_BER = output.BER;
    y1_BER_min(i_snr) = min(nz_BER);        y1_BER_max(i_snr) = max(nz_BER);
    y1_BER_p995(i_snr) = prctile(nz_BER,99.5);
    y1_BER_p99(i_snr) = prctile(nz_BER,99); y1_BER_p95(i_snr) = prctile(nz_BER,95);
    y1_BER_p90(i_snr) = prctile(nz_BER,90); y1_BER_p80(i_snr) = prctile(nz_BER,80);
    y1_BER_p70(i_snr) = prctile(nz_BER,70); y1_BER_p60(i_snr) = prctile(nz_BER,60);
    y1_BER_p50(i_snr) = prctile(nz_BER,50); 
    y1_BER_p40(i_snr) = prctile(nz_BER,40); y1_BER_p30(i_snr) = prctile(nz_BER,30);
    y1_BER_p20(i_snr) = prctile(nz_BER,20); y1_BER_p10(i_snr) = prctile(nz_BER,10);
    y1_BER_p05(i_snr) = prctile(nz_BER, 5); y1_BER_p01(i_snr) = prctile(nz_BER, 1);
    y1_BER_p005(i_snr) = prctile(nz_BER,0.5);
    
    

  end
  yUniform_SINR_min(i_real, :) = mean(y1_SINR_min,1); 
  yUniform_SINR_mean(i_real,:) = mean(y1_SINR_mean,1);
  yUniform_SINR_p995(i_real, :) = mean(y1_SINR_p995,1); 
  yUniform_SINR_p99(i_real, :) = mean(y1_SINR_p99,1); yUniform_SINR_p95(i_real, :) = mean(y1_SINR_p95,1); 
  yUniform_SINR_p90(i_real, :) = mean(y1_SINR_p90,1); yUniform_SINR_p80(i_real, :) = mean(y1_SINR_p80,1);
  yUniform_SINR_p70(i_real, :) = mean(y1_SINR_p70,1); yUniform_SINR_p60(i_real, :) = mean(y1_SINR_p60,1); 
  yUniform_SINR_p50(i_real, :) = mean(y1_SINR_p50,1);
  yUniform_SINR_p40(i_real, :) = mean(y1_SINR_p40,1); yUniform_SINR_p30(i_real, :) = mean(y1_SINR_p30,1);
  yUniform_SINR_p20(i_real, :) = mean(y1_SINR_p20,1); yUniform_SINR_p10(i_real, :) = mean(y1_SINR_p10,1);
  yUniform_SINR_p05(i_real, :) = mean(y1_SINR_p05,1); yUniform_SINR_p01(i_real, :) = mean(y1_SINR_p01,1);
  yUniform_SINR_p005(i_real, :) = mean(y1_SINR_p005,1);

  yUniform_BER_min(i_real, :) = mean(y1_BER_min,1);
  yUniform_BER_max(i_real, :) = mean(y1_BER_max,1);
  yUniform_BER_p995(i_real, :) = mean(y1_BER_p995,1);
  yUniform_BER_p99(i_real, :) = mean(y1_BER_p99,1); yUniform_BER_p95(i_real, :) = mean(y1_BER_p95,1);
  yUniform_BER_p90(i_real, :) = mean(y1_BER_p90,1); yUniform_BER_p80(i_real, :) = mean(y1_BER_p80,1);
  yUniform_BER_p70(i_real, :) = mean(y1_BER_p70,1); yUniform_BER_p60(i_real, :) = mean(y1_BER_p60,1);
  yUniform_BER_p50(i_real, :) = mean(y1_BER_p50,1);
  yUniform_BER_p40(i_real, :) = mean(y1_BER_p40,1); yUniform_BER_p30(i_real, :) = mean(y1_BER_p30,1);
  yUniform_BER_p20(i_real, :) = mean(y1_BER_p20,1); yUniform_BER_p10(i_real, :) = mean(y1_BER_p10,1);
  yUniform_BER_p05(i_real, :) = mean(y1_BER_p05,1); yUniform_BER_p01(i_real, :) = mean(y1_BER_p01,1);
  yUniform_BER_p005(i_real, :) = mean(y1_BER_p005,1);

  % ================ No FH case =====================================
  % X_jam = Jset(:,1:Nh);
  X_jam = Jsets{i_real+1};
  %   Xi_noFH

  oo = transX2Y(Xi_noFH, X_jam, Ns);
  Yi = oo.Yi;
  Y_jam = oo.Y_jam;

  for i_snr = 1:length(SNRdB_list)
    SNRdB = SNRdB_list(i_snr);
    output = simQPSK_Y(Yi, Y_jam, SNRdB);
    nz_SINR = output.SINR_ues;
    
    
    y1_SINR_min(i_snr) = min(nz_SINR);        y1_SINR_mean(i_snr) = mean(nz_SINR);
    y1_SINR_p995(i_snr) = prctile(nz_SINR,99.5); 
    y1_SINR_p99(i_snr) = prctile(nz_SINR,99); y1_SINR_p95(i_snr) = prctile(nz_SINR,95);
    y1_SINR_p90(i_snr) = prctile(nz_SINR,90); y1_SINR_p80(i_snr) = prctile(nz_SINR,80);
    y1_SINR_p70(i_snr) = prctile(nz_SINR,70); y1_SINR_p60(i_snr) = prctile(nz_SINR,60);
    y1_SINR_p50(i_snr) = prctile(nz_SINR,50); 
    y1_SINR_p40(i_snr) = prctile(nz_SINR,40); y1_SINR_p30(i_snr) = prctile(nz_SINR,30);
    y1_SINR_p20(i_snr) = prctile(nz_SINR,20); y1_SINR_p10(i_snr) = prctile(nz_SINR,10);
    y1_SINR_p05(i_snr) = prctile(nz_SINR, 5); y1_SINR_p01(i_snr) = prctile(nz_SINR, 1); 
    y1_SINR_p005(i_snr) = prctile(nz_SINR,0.5);

    output = simBER_givenYiYjam(Yi, Y_jam, SNRdB);
    nz_BER = output.BER;
    y1_BER_min(i_snr) = min(nz_BER);        y1_BER_max(i_snr) = max(nz_BER);
    y1_BER_p995(i_snr) = prctile(nz_BER,99.5);
    y1_BER_p99(i_snr) = prctile(nz_BER,99); y1_BER_p95(i_snr) = prctile(nz_BER,95);
    y1_BER_p90(i_snr) = prctile(nz_BER,90); y1_BER_p80(i_snr) = prctile(nz_BER,80);
    y1_BER_p70(i_snr) = prctile(nz_BER,70); y1_BER_p60(i_snr) = prctile(nz_BER,60);
    y1_BER_p50(i_snr) = prctile(nz_BER,50); 
    y1_BER_p40(i_snr) = prctile(nz_BER,40); y1_BER_p30(i_snr) = prctile(nz_BER,30);
    y1_BER_p20(i_snr) = prctile(nz_BER,20); y1_BER_p10(i_snr) = prctile(nz_BER,10);
    y1_BER_p05(i_snr) = prctile(nz_BER, 5); y1_BER_p01(i_snr) = prctile(nz_BER, 1);
    y1_BER_p005(i_snr) = prctile(nz_BER,0.5);
    
    

  end
  yNoFH_SINR_min(i_real, :) = mean(y1_SINR_min,1); 
  yNoFH_SINR_mean(i_real,:) = mean(y1_SINR_mean,1);
  yNoFH_SINR_p995(i_real, :) = mean(y1_SINR_p995,1); 
  yNoFH_SINR_p99(i_real, :) = mean(y1_SINR_p99,1); yNoFH_SINR_p95(i_real, :) = mean(y1_SINR_p95,1); 
  yNoFH_SINR_p90(i_real, :) = mean(y1_SINR_p90,1); yNoFH_SINR_p80(i_real, :) = mean(y1_SINR_p80,1);
  yNoFH_SINR_p70(i_real, :) = mean(y1_SINR_p70,1); yNoFH_SINR_p60(i_real, :) = mean(y1_SINR_p60,1); 
  yNoFH_SINR_p50(i_real, :) = mean(y1_SINR_p50,1);
  yNoFH_SINR_p40(i_real, :) = mean(y1_SINR_p40,1); yNoFH_SINR_p30(i_real, :) = mean(y1_SINR_p30,1);
  yNoFH_SINR_p20(i_real, :) = mean(y1_SINR_p20,1); yNoFH_SINR_p10(i_real, :) = mean(y1_SINR_p10,1);
  yNoFH_SINR_p05(i_real, :) = mean(y1_SINR_p05,1); yNoFH_SINR_p01(i_real, :) = mean(y1_SINR_p01,1);
  yNoFH_SINR_p005(i_real, :) = mean(y1_SINR_p005,1);

  yNoFH_BER_min(i_real, :) = mean(y1_BER_min,1);
  yNoFH_BER_max(i_real, :) = mean(y1_BER_max,1);
  yNoFH_BER_p995(i_real, :) = mean(y1_BER_p995,1);
  yNoFH_BER_p99(i_real, :) = mean(y1_BER_p99,1); yNoFH_BER_p95(i_real, :) = mean(y1_BER_p95,1);
  yNoFH_BER_p90(i_real, :) = mean(y1_BER_p90,1); yNoFH_BER_p80(i_real, :) = mean(y1_BER_p80,1);
  yNoFH_BER_p70(i_real, :) = mean(y1_BER_p70,1); yNoFH_BER_p60(i_real, :) = mean(y1_BER_p60,1);
  yNoFH_BER_p50(i_real, :) = mean(y1_BER_p50,1);
  yNoFH_BER_p40(i_real, :) = mean(y1_BER_p40,1); yNoFH_BER_p30(i_real, :) = mean(y1_BER_p30,1);
  yNoFH_BER_p20(i_real, :) = mean(y1_BER_p20,1); yNoFH_BER_p10(i_real, :) = mean(y1_BER_p10,1);
  yNoFH_BER_p05(i_real, :) = mean(y1_BER_p05,1); yNoFH_BER_p01(i_real, :) = mean(y1_BER_p01,1);
  yNoFH_BER_p005(i_real, :) = mean(y1_BER_p005,1);

  
  
  % ================ No Jamming case =====================================
  % X_jam = Jset(:,1:Nh);
  X_jam = Jsets{i_real+1};
  % oo = transX2Y(Xi_Stat, X_jam, Ns);
  oo = transX2Y(Xi_Det, X_jam, Ns);
  Yi = oo.Yi;
  Y_jam = zeros(Ns, Nh); % oo.Y_jam;

%   SNRdB_list = [-2 0 2 4 6 7];
  for i_snr = 1:length(SNRdB_list)
    SNRdB = SNRdB_list(i_snr);
    output = simQPSK_Y(Yi, Y_jam, SNRdB);
    nz_SINR = output.SINR_ues;
    
    
    y1_SINR_min(i_snr) = min(nz_SINR);        y1_SINR_mean(i_snr) = mean(nz_SINR);
    y1_SINR_p995(i_snr) = prctile(nz_SINR,99.5); 
    y1_SINR_p99(i_snr) = prctile(nz_SINR,99); y1_SINR_p95(i_snr) = prctile(nz_SINR,95);
    y1_SINR_p90(i_snr) = prctile(nz_SINR,90); y1_SINR_p80(i_snr) = prctile(nz_SINR,80);
    y1_SINR_p70(i_snr) = prctile(nz_SINR,70); y1_SINR_p60(i_snr) = prctile(nz_SINR,60);
    y1_SINR_p50(i_snr) = prctile(nz_SINR,50); 
    y1_SINR_p40(i_snr) = prctile(nz_SINR,40); y1_SINR_p30(i_snr) = prctile(nz_SINR,30);
    y1_SINR_p20(i_snr) = prctile(nz_SINR,20); y1_SINR_p10(i_snr) = prctile(nz_SINR,10);
    y1_SINR_p05(i_snr) = prctile(nz_SINR, 5); y1_SINR_p01(i_snr) = prctile(nz_SINR, 1); 
    y1_SINR_p005(i_snr) = prctile(nz_SINR,0.5);
    
  
    output = simBER_givenYiYjam(Yi, Y_jam, SNRdB);
    nz_BER = output.BER;
    y1_BER_min(i_snr) = min(nz_BER);        y1_BER_max(i_snr) = max(nz_BER);
    y1_BER_p995(i_snr) = prctile(nz_BER,99.5);
    y1_BER_p99(i_snr) = prctile(nz_BER,99); y1_BER_p95(i_snr) = prctile(nz_BER,95);
    y1_BER_p90(i_snr) = prctile(nz_BER,90); y1_BER_p80(i_snr) = prctile(nz_BER,80);
    y1_BER_p70(i_snr) = prctile(nz_BER,70); y1_BER_p60(i_snr) = prctile(nz_BER,60);
    y1_BER_p50(i_snr) = prctile(nz_BER,50); 
    y1_BER_p40(i_snr) = prctile(nz_BER,40); y1_BER_p30(i_snr) = prctile(nz_BER,30);
    y1_BER_p20(i_snr) = prctile(nz_BER,20); y1_BER_p10(i_snr) = prctile(nz_BER,10);
    y1_BER_p05(i_snr) = prctile(nz_BER, 5); y1_BER_p01(i_snr) = prctile(nz_BER, 1);
    y1_BER_p005(i_snr) = prctile(nz_BER,0.5);
    

  end
  yNoJam_SINR_min(i_real, :) = mean(y1_SINR_min,1); 
  yNoJam_SINR_mean(i_real,:) = mean(y1_SINR_mean,1);
  yNoJam_SINR_p995(i_real, :) = mean(y1_SINR_p995,1); 
  yNoJam_SINR_p99(i_real, :) = mean(y1_SINR_p99,1); yNoJam_SINR_p95(i_real, :) = mean(y1_SINR_p95,1); 
  yNoJam_SINR_p90(i_real, :) = mean(y1_SINR_p90,1); yNoJam_SINR_p80(i_real, :) = mean(y1_SINR_p80,1);
  yNoJam_SINR_p70(i_real, :) = mean(y1_SINR_p70,1); yNoJam_SINR_p60(i_real, :) = mean(y1_SINR_p60,1); 
  yNoJam_SINR_p50(i_real, :) = mean(y1_SINR_p50,1);
  yNoJam_SINR_p40(i_real, :) = mean(y1_SINR_p40,1); yNoJam_SINR_p30(i_real, :) = mean(y1_SINR_p30,1);
  yNoJam_SINR_p20(i_real, :) = mean(y1_SINR_p20,1); yNoJam_SINR_p10(i_real, :) = mean(y1_SINR_p10,1);
  yNoJam_SINR_p05(i_real, :) = mean(y1_SINR_p05,1); yNoJam_SINR_p01(i_real, :) = mean(y1_SINR_p01,1);
  yNoJam_SINR_p005(i_real, :) = mean(y1_SINR_p005,1);
  
  
  yNoJam_BER_min(i_real, :) = mean(y1_BER_min,1);
  yNoJam_BER_max(i_real, :) = mean(y1_BER_max,1);
  yNoJam_BER_p995(i_real, :) = mean(y1_BER_p995,1);
  yNoJam_BER_p99(i_real, :) = mean(y1_BER_p99,1); yNoJam_BER_p95(i_real, :) = mean(y1_BER_p95,1);
  yNoJam_BER_p90(i_real, :) = mean(y1_BER_p90,1); yNoJam_BER_p80(i_real, :) = mean(y1_BER_p80,1);
  yNoJam_BER_p70(i_real, :) = mean(y1_BER_p70,1); yNoJam_BER_p60(i_real, :) = mean(y1_BER_p60,1);
  yNoJam_BER_p50(i_real, :) = mean(y1_BER_p50,1);
  yNoJam_BER_p40(i_real, :) = mean(y1_BER_p40,1); yNoJam_BER_p30(i_real, :) = mean(y1_BER_p30,1);
  yNoJam_BER_p20(i_real, :) = mean(y1_BER_p20,1); yNoJam_BER_p10(i_real, :) = mean(y1_BER_p10,1);
  yNoJam_BER_p05(i_real, :) = mean(y1_BER_p05,1); yNoJam_BER_p01(i_real, :) = mean(y1_BER_p01,1);
  yNoJam_BER_p005(i_real, :) = mean(y1_BER_p005,1);
  
end


% Average over all realizetions
yDet_noncausal_SINR_min = mean(yDet_noncausal_SINR_min,1);
yDet_noncausal_SINR_mean = mean(yDet_noncausal_SINR_mean,1);
yDet_noncausal_SINR_p995 = mean(yDet_noncausal_SINR_p995,1);
yDet_noncausal_SINR_p99 = mean(yDet_noncausal_SINR_p99,1);
yDet_noncausal_SINR_p95 = mean(yDet_noncausal_SINR_p95,1);
yDet_noncausal_SINR_p90 = mean(yDet_noncausal_SINR_p90,1);
yDet_noncausal_SINR_p80 = mean(yDet_noncausal_SINR_p80,1);
yDet_noncausal_SINR_p70 = mean(yDet_noncausal_SINR_p70,1);
yDet_noncausal_SINR_p60 = mean(yDet_noncausal_SINR_p60,1);
yDet_noncausal_SINR_p50 = mean(yDet_noncausal_SINR_p50,1);
yDet_noncausal_SINR_p40 = mean(yDet_noncausal_SINR_p40,1);
yDet_noncausal_SINR_p30 = mean(yDet_noncausal_SINR_p30,1);
yDet_noncausal_SINR_p20 = mean(yDet_noncausal_SINR_p20,1);
yDet_noncausal_SINR_p10 = mean(yDet_noncausal_SINR_p10,1);
yDet_noncausal_SINR_p05 = mean(yDet_noncausal_SINR_p05,1);
yDet_noncausal_SINR_p01 = mean(yDet_noncausal_SINR_p01,1);
yDet_noncausal_SINR_p005 = mean(yDet_noncausal_SINR_p005,1);

yDet_noncausal_BER_min = mean(yDet_noncausal_BER_min,1);
yDet_noncausal_BER_p995 = mean(yDet_noncausal_BER_p995,1);
yDet_noncausal_BER_p99 = mean(yDet_noncausal_BER_p99,1);
yDet_noncausal_BER_p95 = mean(yDet_noncausal_BER_p95,1);
yDet_noncausal_BER_p90 = mean(yDet_noncausal_BER_p90,1);
yDet_noncausal_BER_p80 = mean(yDet_noncausal_BER_p80,1);
yDet_noncausal_BER_p70 = mean(yDet_noncausal_BER_p70,1);
yDet_noncausal_BER_p60 = mean(yDet_noncausal_BER_p60,1);
yDet_noncausal_BER_p50 = mean(yDet_noncausal_BER_p50,1);
yDet_noncausal_BER_p40 = mean(yDet_noncausal_BER_p40,1);
yDet_noncausal_BER_p30 = mean(yDet_noncausal_BER_p30,1);
yDet_noncausal_BER_p20 = mean(yDet_noncausal_BER_p20,1);
yDet_noncausal_BER_p10 = mean(yDet_noncausal_BER_p10,1);
yDet_noncausal_BER_p05 = mean(yDet_noncausal_BER_p05,1);
yDet_noncausal_BER_p01 = mean(yDet_noncausal_BER_p01,1);
yDet_noncausal_BER_p005 = mean(yDet_noncausal_BER_p005,1);
yDet_noncausal_BER_max = mean(yDet_noncausal_BER_max,1);

yDet_causal_SINR_min = mean(yDet_causal_SINR_min,1);
yDet_causal_SINR_mean = mean(yDet_causal_SINR_mean,1);
yDet_causal_SINR_p995 = mean(yDet_causal_SINR_p995,1);
yDet_causal_SINR_p99 = mean(yDet_causal_SINR_p99,1);
yDet_causal_SINR_p95 = mean(yDet_causal_SINR_p95,1);
yDet_causal_SINR_p90 = mean(yDet_causal_SINR_p90,1);
yDet_causal_SINR_p80 = mean(yDet_causal_SINR_p80,1);
yDet_causal_SINR_p70 = mean(yDet_causal_SINR_p70,1);
yDet_causal_SINR_p60 = mean(yDet_causal_SINR_p60,1);
yDet_causal_SINR_p50 = mean(yDet_causal_SINR_p50,1);
yDet_causal_SINR_p40 = mean(yDet_causal_SINR_p40,1);
yDet_causal_SINR_p30 = mean(yDet_causal_SINR_p30,1);
yDet_causal_SINR_p20 = mean(yDet_causal_SINR_p20,1);
yDet_causal_SINR_p10 = mean(yDet_causal_SINR_p10,1);
yDet_causal_SINR_p05 = mean(yDet_causal_SINR_p05,1);
yDet_causal_SINR_p01 = mean(yDet_causal_SINR_p01,1);
yDet_causal_SINR_p005 = mean(yDet_causal_SINR_p005,1);

yDet_causal_BER_min= mean(yDet_causal_BER_min,1);
yDet_causal_BER_p995 = mean(yDet_causal_BER_p995,1);
yDet_causal_BER_p99 = mean(yDet_causal_BER_p99,1);
yDet_causal_BER_p95 = mean(yDet_causal_BER_p95,1);
yDet_causal_BER_p90 = mean(yDet_causal_BER_p90,1);
yDet_causal_BER_p80 = mean(yDet_causal_BER_p80,1);
yDet_causal_BER_p70 = mean(yDet_causal_BER_p70,1);
yDet_causal_BER_p60 = mean(yDet_causal_BER_p60,1);
yDet_causal_BER_p50 = mean(yDet_causal_BER_p50,1);
yDet_causal_BER_p40 = mean(yDet_causal_BER_p40,1);
yDet_causal_BER_p30 = mean(yDet_causal_BER_p30,1);
yDet_causal_BER_p20 = mean(yDet_causal_BER_p20,1);
yDet_causal_BER_p10 = mean(yDet_causal_BER_p10,1);
yDet_causal_BER_p05 = mean(yDet_causal_BER_p05,1);
yDet_causal_BER_p01 = mean(yDet_causal_BER_p01,1);
yDet_causal_BER_p005 = mean(yDet_causal_BER_p005,1);
yDet_causal_BER_max = mean(yDet_causal_BER_max,1);

yStat_SINR_min = mean(yStat_SINR_min,1);
yStat_SINR_mean= mean(yStat_SINR_mean,1);
yStat_SINR_p995 = mean(yStat_SINR_p995,1);
yStat_SINR_p99 = mean(yStat_SINR_p99,1);
yStat_SINR_p95 = mean(yStat_SINR_p95,1);
yStat_SINR_p90 = mean(yStat_SINR_p90,1);
yStat_SINR_p80 = mean(yStat_SINR_p80,1);
yStat_SINR_p70 = mean(yStat_SINR_p70,1);
yStat_SINR_p60 = mean(yStat_SINR_p60,1);
yStat_SINR_p50 = mean(yStat_SINR_p50,1);
yStat_SINR_p40 = mean(yStat_SINR_p40,1);
yStat_SINR_p30 = mean(yStat_SINR_p30,1);
yStat_SINR_p20 = mean(yStat_SINR_p20,1);
yStat_SINR_p10 = mean(yStat_SINR_p10,1);
yStat_SINR_p05 = mean(yStat_SINR_p05,1);
yStat_SINR_p01 = mean(yStat_SINR_p01,1);
yStat_SINR_p005 = mean(yStat_SINR_p005,1);

yStat_BER_min = mean(yStat_BER_min,1);
yStat_BER_p995 = mean(yStat_BER_p995,1);
yStat_BER_p99 = mean(yStat_BER_p99,1);
yStat_BER_p95 = mean(yStat_BER_p95,1);
yStat_BER_p90 = mean(yStat_BER_p90,1);
yStat_BER_p80 = mean(yStat_BER_p80,1);
yStat_BER_p70 = mean(yStat_BER_p70,1);
yStat_BER_p60 = mean(yStat_BER_p60,1);
yStat_BER_p50 = mean(yStat_BER_p50,1);
yStat_BER_p40 = mean(yStat_BER_p40,1);
yStat_BER_p30 = mean(yStat_BER_p30,1);
yStat_BER_p20 = mean(yStat_BER_p20,1);
yStat_BER_p10 = mean(yStat_BER_p10,1);
yStat_BER_p05 = mean(yStat_BER_p05,1);
yStat_BER_p01 = mean(yStat_BER_p01,1);
yStat_BER_p005 = mean(yStat_BER_p005,1);
yStat_BER_max = mean(yStat_BER_max,1);

yStat_m2_SINR_min = mean(yStat_m2_SINR_min,1);
yStat_m2_SINR_mean= mean(yStat_m2_SINR_mean,1);
yStat_m2_SINR_p995 = mean(yStat_m2_SINR_p995,1);
yStat_m2_SINR_p99 = mean(yStat_m2_SINR_p99,1);
yStat_m2_SINR_p95 = mean(yStat_m2_SINR_p95,1);
yStat_m2_SINR_p90 = mean(yStat_m2_SINR_p90,1);
yStat_m2_SINR_p80 = mean(yStat_m2_SINR_p80,1);
yStat_m2_SINR_p70 = mean(yStat_m2_SINR_p70,1);
yStat_m2_SINR_p60 = mean(yStat_m2_SINR_p60,1);
yStat_m2_SINR_p50 = mean(yStat_m2_SINR_p50,1);
yStat_m2_SINR_p40 = mean(yStat_m2_SINR_p40,1);
yStat_m2_SINR_p30 = mean(yStat_m2_SINR_p30,1);
yStat_m2_SINR_p20 = mean(yStat_m2_SINR_p20,1);
yStat_m2_SINR_p10 = mean(yStat_m2_SINR_p10,1);
yStat_m2_SINR_p05 = mean(yStat_m2_SINR_p05,1);
yStat_m2_SINR_p01 = mean(yStat_m2_SINR_p01,1);
yStat_m2_SINR_p005 = mean(yStat_m2_SINR_p005,1);

yStat_m2_BER_min = mean(yStat_m2_BER_min,1);
yStat_m2_BER_p995 = mean(yStat_m2_BER_p995,1);
yStat_m2_BER_p99 = mean(yStat_m2_BER_p99,1);
yStat_m2_BER_p95 = mean(yStat_m2_BER_p95,1);
yStat_m2_BER_p90 = mean(yStat_m2_BER_p90,1);
yStat_m2_BER_p80 = mean(yStat_m2_BER_p80,1);
yStat_m2_BER_p70 = mean(yStat_m2_BER_p70,1);
yStat_m2_BER_p60 = mean(yStat_m2_BER_p60,1);
yStat_m2_BER_p50 = mean(yStat_m2_BER_p50,1);
yStat_m2_BER_p40 = mean(yStat_m2_BER_p40,1);
yStat_m2_BER_p30 = mean(yStat_m2_BER_p30,1);
yStat_m2_BER_p20 = mean(yStat_m2_BER_p20,1);
yStat_m2_BER_p10 = mean(yStat_m2_BER_p10,1);
yStat_m2_BER_p05 = mean(yStat_m2_BER_p05,1);
yStat_m2_BER_p01 = mean(yStat_m2_BER_p01,1);
yStat_m2_BER_p005 = mean(yStat_m2_BER_p005,1);
yStat_m2_BER_max = mean(yStat_m2_BER_max,1);

yUniform_SINR_min = mean(yUniform_SINR_min,1);
yUniform_SINR_mean= mean(yUniform_SINR_mean,1);
yUniform_SINR_p995 = mean(yUniform_SINR_p995,1);
yUniform_SINR_p99 = mean(yUniform_SINR_p99,1);
yUniform_SINR_p95 = mean(yUniform_SINR_p95,1);
yUniform_SINR_p90 = mean(yUniform_SINR_p90,1);
yUniform_SINR_p80 = mean(yUniform_SINR_p80,1);
yUniform_SINR_p70 = mean(yUniform_SINR_p70,1);
yUniform_SINR_p60 = mean(yUniform_SINR_p60,1);
yUniform_SINR_p50 = mean(yUniform_SINR_p50,1);
yUniform_SINR_p40 = mean(yUniform_SINR_p40,1);
yUniform_SINR_p30 = mean(yUniform_SINR_p30,1);
yUniform_SINR_p20 = mean(yUniform_SINR_p20,1);
yUniform_SINR_p10 = mean(yUniform_SINR_p10,1);
yUniform_SINR_p05 = mean(yUniform_SINR_p05,1);
yUniform_SINR_p01 = mean(yUniform_SINR_p01,1);
yUniform_SINR_p005 = mean(yUniform_SINR_p005,1);

yUniform_BER_min = mean(yUniform_BER_min,1);
yUniform_BER_p995 = mean(yUniform_BER_p995,1);
yUniform_BER_p99 = mean(yUniform_BER_p99,1);
yUniform_BER_p95 = mean(yUniform_BER_p95,1);
yUniform_BER_p90 = mean(yUniform_BER_p90,1);
yUniform_BER_p80 = mean(yUniform_BER_p80,1);
yUniform_BER_p70 = mean(yUniform_BER_p70,1);
yUniform_BER_p60 = mean(yUniform_BER_p60,1);
yUniform_BER_p50 = mean(yUniform_BER_p50,1);
yUniform_BER_p40 = mean(yUniform_BER_p40,1);
yUniform_BER_p30 = mean(yUniform_BER_p30,1);
yUniform_BER_p20 = mean(yUniform_BER_p20,1);
yUniform_BER_p10 = mean(yUniform_BER_p10,1);
yUniform_BER_p05 = mean(yUniform_BER_p05,1);
yUniform_BER_p01 = mean(yUniform_BER_p01,1);
yUniform_BER_p005 = mean(yUniform_BER_p005,1);
yUniform_BER_max = mean(yUniform_BER_max,1);



yNoFH_SINR_min = mean(yNoFH_SINR_min,1);
yNoFH_SINR_mean= mean(yNoFH_SINR_mean,1);
yNoFH_SINR_p995 = mean(yNoFH_SINR_p995,1);
yNoFH_SINR_p99 = mean(yNoFH_SINR_p99,1);
yNoFH_SINR_p95 = mean(yNoFH_SINR_p95,1);
yNoFH_SINR_p90 = mean(yNoFH_SINR_p90,1);
yNoFH_SINR_p80 = mean(yNoFH_SINR_p80,1);
yNoFH_SINR_p70 = mean(yNoFH_SINR_p70,1);
yNoFH_SINR_p60 = mean(yNoFH_SINR_p60,1);
yNoFH_SINR_p50 = mean(yNoFH_SINR_p50,1);
yNoFH_SINR_p40 = mean(yNoFH_SINR_p40,1);
yNoFH_SINR_p30 = mean(yNoFH_SINR_p30,1);
yNoFH_SINR_p20 = mean(yNoFH_SINR_p20,1);
yNoFH_SINR_p10 = mean(yNoFH_SINR_p10,1);
yNoFH_SINR_p05 = mean(yNoFH_SINR_p05,1);
yNoFH_SINR_p01 = mean(yNoFH_SINR_p01,1);
yNoFH_SINR_p005 = mean(yNoFH_SINR_p005,1);

yNoFH_BER_min = mean(yUniform_BER_min,1);
yNoFH_BER_p995 = mean(yNoFH_BER_p995,1);
yNoFH_BER_p99 = mean(yNoFH_BER_p99,1);
yNoFH_BER_p95 = mean(yNoFH_BER_p95,1);
yNoFH_BER_p90 = mean(yNoFH_BER_p90,1);
yNoFH_BER_p80 = mean(yNoFH_BER_p80,1);
yNoFH_BER_p70 = mean(yNoFH_BER_p70,1);
yNoFH_BER_p60 = mean(yNoFH_BER_p60,1);
yNoFH_BER_p50 = mean(yNoFH_BER_p50,1);
yNoFH_BER_p40 = mean(yNoFH_BER_p40,1);
yNoFH_BER_p30 = mean(yNoFH_BER_p30,1);
yNoFH_BER_p20 = mean(yNoFH_BER_p20,1);
yNoFH_BER_p10 = mean(yNoFH_BER_p10,1);
yNoFH_BER_p05 = mean(yNoFH_BER_p05,1);
yNoFH_BER_p01 = mean(yNoFH_BER_p01,1);
yNoFH_BER_p005 = mean(yNoFH_BER_p005,1);
yNoFH_BER_max = mean(yNoFH_BER_max,1);


yNoJam_SINR_min = mean(yNoJam_SINR_min,1);
yNoJam_SINR_mean= mean(yNoJam_SINR_mean,1);
yNoJam_SINR_p995 = mean(yNoJam_SINR_p995,1);
yNoJam_SINR_p99 = mean(yNoJam_SINR_p99,1);
yNoJam_SINR_p95 = mean(yNoJam_SINR_p95,1);
yNoJam_SINR_p90 = mean(yNoJam_SINR_p90,1);
yNoJam_SINR_p80 = mean(yNoJam_SINR_p80,1);
yNoJam_SINR_p70 = mean(yNoJam_SINR_p70,1);
yNoJam_SINR_p60 = mean(yNoJam_SINR_p60,1);
yNoJam_SINR_p50 = mean(yNoJam_SINR_p50,1);
yNoJam_SINR_p40 = mean(yNoJam_SINR_p40,1);
yNoJam_SINR_p30 = mean(yNoJam_SINR_p30,1);
yNoJam_SINR_p20 = mean(yNoJam_SINR_p20,1);
yNoJam_SINR_p10 = mean(yNoJam_SINR_p10,1);
yNoJam_SINR_p05 = mean(yNoJam_SINR_p05,1);
yNoJam_SINR_p01 = mean(yNoJam_SINR_p01,1);
yNoJam_SINR_p005 = mean(yNoJam_SINR_p005,1);

yNoJam_BER_min = mean(yNoJam_BER_min,1);
yNoJam_BER_p995 = mean(yNoJam_BER_p995,1);
yNoJam_BER_p99 = mean(yNoJam_BER_p99,1);
yNoJam_BER_p95 = mean(yNoJam_BER_p95,1);
yNoJam_BER_p90 = mean(yNoJam_BER_p90,1);
yNoJam_BER_p80 = mean(yNoJam_BER_p80,1);
yNoJam_BER_p70 = mean(yNoJam_BER_p70,1);
yNoJam_BER_p60 = mean(yNoJam_BER_p60,1);
yNoJam_BER_p50 = mean(yNoJam_BER_p50,1);
yNoJam_BER_p40 = mean(yNoJam_BER_p40,1);
yNoJam_BER_p30 = mean(yNoJam_BER_p30,1);
yNoJam_BER_p20 = mean(yNoJam_BER_p20,1);
yNoJam_BER_p10 = mean(yNoJam_BER_p10,1);
yNoJam_BER_p05 = mean(yNoJam_BER_p05,1);
yNoJam_BER_p01 = mean(yNoJam_BER_p01,1);
yNoJam_BER_p005 = mean(yNoJam_BER_p005,1);
yNoJam_BER_max = mean(yNoJam_BER_max,1);



fprintf("------------- End Data Generation --------------------------\n")
save fh_low2SNR.mat
% save fh_mid4SNR.mat
% save fh_highSNR.mat
% save fh_testSNR.mat

genFigures;





















%===========================================================================
function [output] = dist_H(X1, X2)
  output = length( find( X1-X2 ~= 0 ) );
end


function [output] = genPattern(pID, Ns, Nu)
  % Given one particular permutation index,
  % Generate one particular pattern for all UEs
    for uID = 1:Nu
      pattern_u = mod(floor( (pID-1)/ (Ns^(uID-1))), Ns) + 1;
      X_set(uID,:) = pattern_u;
    end
    output = X_set;
end




function [output] = transX2Y(Xi, X_jam, Ns)
  Nu = size(Xi,1);
  Nh = size(Xi,2);
  Nj = size(X_jam,1);

  for i=1:Nu
    Y_tmp = zeros(Ns, Nh);
    anz = find(Xi(i,:)>0);
    for l = 1:length(anz)
      Y_tmp( Xi(i,anz(l)), anz(l) ) = 1;
    end
    Yi{i} = Y_tmp;
  end

  Y_jam = zeros(Ns,Nh);
  for n=1:Nh
    for l=1:Nj
      fj = X_jam(l,n);
      if (fj > 0)
        Y_jam(fj,n) = 1;
      end
    end
  end

  output.Yi = Yi;
  output.Y_jam = Y_jam;

end

function [output] = simQPSK_Y(Yi, Y_jam, SNRdB)
  SNRlin = 10.^(SNRdB/10);
  Nu = length(Yi);
  Yitmp = Yi{1};
  Ns = size(Yitmp,1);
  Nh = size(Yitmp,2);

  Y_total = Y_jam;
  for i=1:Nu
    Y_total = Y_total + Yi{i};
  end

  for i=1:Nu
    sumSINR = 0;
    count = 0;
    for n=1:Nh
      yi_mat = Yi{i};
      yi_n = yi_mat(:,n);
      y_total = Y_total(:,n);

      Anz = find(yi_n > 0);
      Dnz = find( (y_total-yi_n).*yi_n >0);
      if (length(Anz)>0)
        nhit = length(Dnz) / length(Anz);
        sumSINR = sumSINR + 10*log10(SNRlin / (1+ SNRlin*nhit));
        count = count + 1;
      end
    end
    SINR_ues(i) = sumSINR / count;
  end

  output.SINR_ues =SINR_ues; 
end

function [pj] = genTrue_pj(pJ_true)
  mu = pJ_true.mu;
  sigma = pJ_true.sigma;
  a = pJ_true.a;
  b = pJ_true.b;
  % n = pJ_true.n;

  deno= fun_Phi( (b-mu)/ sigma ) - fun_Phi( (a-mu)/sigma );
  for i=a:b
      
    pj(i) = (1/sigma) * (1/sqrt(2*pi)) * exp(-0.5*( (i-mu)/sigma )^2) / (deno);
  end
end

function [pmf] = genPMF_pj( samples, max_b )
  n = size(samples,1);
  
  pmf = accumarray( samples(:),1)./numel(samples);
  pmf = pmf./sum(pmf);
  if (size(pmf,1) < max_b)
    pmf(max_b,1) = 0;
  end
end

function [output] = genJamPattern(pJ_true)
  mu = pJ_true.mu;
  sigma = pJ_true.sigma;
  a = pJ_true.a;
  b = pJ_true.b;
  n = pJ_true.n;

  U = rand(n,1);
  tmp = U * ( fun_Phi( (b-mu)/ sigma) - fun_Phi( (a-mu)/sigma ) ) + fun_Phi( (a-mu)/sigma );
  output = sqrt(2) * erfinv(-1 + 2 * tmp) * sigma + mu;
  output = round(output);
end

function [output] = fun_Phi(x)
  output = 0.5 * (1 + erf(x / sqrt(2)) );
end


function [output] = simBER_givenYiYjam(Yi, Yjam, snrdB)
  % Simulate BER performance for whole Nh slot

  Ns = size(Yjam,1);
  Nh = size(Yjam,2);
  Nu = length(Yi);
  Nsym = 1e2; % Number of symbols transmitted per FH interval 
  
  Y_supp = zeros(Ns, Nh); % Support matrix to indicate UE transmission

  for n=1:Nh
    
    b1 = rand(Nu+1, Nsym) > 0.5;
    b2 = rand(Nu+1, Nsym) > 0.5;
  
    % QPSK sumbol mapping
    I = (2*b1) - 1;
    Q = (2*b2) - 1;
    S = I + 1j*Q;

    SNRlin = 10.^(snrdB/10);
    N0 = 1./SNRlin;   % Noise Variance

    noise = sqrt(N0/2) * (randn(Ns,Nsym) + 1j*randn(Ns,Nsym)); % AWGN noise (row: subcarrier, col: symbol)
    sig_Rx = zeros(Ns, Nsym) + noise;

    for ueI = 1:Nu
      % each user, extract its subcarrier index
      Yi_tmp = Yi{ueI};
      yi = Yi_tmp(:,n);

      setK = find(yi > 0);
      for k_idx = 1:length(setK)
        k = setK(k_idx);

        % superpost the QPSK signals into this subcarrier
        sig_Rx(k,:) = sig_Rx(k,:) + S(ueI,:);
      end
    end

    % Add the jammer's signals
    y_jam = Yjam(:,n);
    setK = find(y_jam > 0);
    for k_idx = 1:length(setK)
        k = setK(k_idx);
        sig_Rx(k,:) = sig_Rx(k,:) + S(Nu+1,:);
    end

    % Decode the symbol for each UE
    % BER(ueI,n) = 0;
    BER_array = [];
    for ueI = 1:Nu
        Yi_tmp = Yi{ueI};
        yi = Yi_tmp(:,n);

        setK = find(yi > 0);
        ber_oneUE = 0;
        for k_idx = 1:length(setK)
            k = setK(k_idx);

            sig_I = real(sig_Rx(k,:)); % I component
            sig_Q = imag(sig_Rx(k,:)); % Q component
        
            bld_I = sig_I > 0; % I decision 
            bld_Q = sig_Q > 0; % Q decision
        
            b1_error = (bld_I ~= b1(ueI,:)); % Inphase bit error
            b2_error = (bld_Q ~= b2(ueI,:)); % Quadrature bit error
        
            Error_bit = sum(b1_error(:)) + sum(b2_error(:)); % Total bit error
            % BER(ueI,n) = sum(Error_bit)/(2*Nsym); % Simulated BER
            % BER(ueI,n) = BER(ueI,n) + sum(Error_bit)/1; % Simulated BER
            ber_oneUE = ber_oneUE + sum(Error_bit);

        end
        if (length(setK) > 0)
          % BER_array = [BER_array, ber_oneUE/ (2*Nsym*length(setK))];
          BER(ueI, n) =  ber_oneUE/ (2*Nsym*length(setK)); 
          Y_supp(ueI, n) = 1;
        else
            BER(ueI, n) = 0;
            Y_supp(ueI, n) = 0;
        end
    end
    % BER(:,n) = sort(BER_array);
  end
  %   output.BER = mean(BER,2);
  output.BER = sum(BER,2)./sum(Y_supp,2);
end


function [output] = genX_LinearFH(Ns,Nh,Nu)
    
    Xi = zeros(Ns, Nh);
    nextUE = 1;

    for n=1:Nh
        for i=1:Nu
            ii = mod(nextUE - 1 + i-1, Nu) + 1;

            if(i <= Ns)
                Xi(ii,n) = i;
                nextUEtmp = ii+1;
            else
                Xi(ii,n) = 0;
            end
        end
        nextUE = nextUEtmp;
    end
    output= Xi;
end

function [output] = genX_noFH(Ns,Nh,Nu)
    
    Xi = zeros(Ns, Nh);
    nextUE = 1;

    for n=1:2:Nh
        for i=1:Ns
            Xi(i,n) = i;
        end
    end
    for n=2:2:Nh
        for i=1:Ns
            Xi(min(Ns+i, Nu),n) = i;
        end
    end
    output= Xi;
end