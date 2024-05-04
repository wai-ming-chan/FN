clc; clear all; close all;
% load fh.mat
rng(12345)

Ns = 1620; % total number of sc
Nh = 14; % number of OFDM symbol in one slot time
Nu1 = 32; % number of URLLC UE
Nu2 = 3000; % number of mMTC UE
Nu = Nu1 + Nu2;
Nj = 32; % number of jammed sc
Nj_band = 512; % number of sc in jamming band
Ni1 = 1; % demanded number of sc for URLLC UE
Ni2 = 1; % demanded number of sc for mMTC UE

for i=1:Nu
  Yi{i} = zeros(Ns, Nh);
end

% Deterministic Jamming FH pattern case
J = zeros(Ns, Nh);
mid_jb = ceil(Ns/2); 
Nj_lb = mid_jb - Nj_band /2;
Nj_ub = mid_jb + Nj_band /2;

for n=1:Nh
  jammed_sc_idx = Nj_lb + randperm(Nj_band, Nj);
  J(jammed_sc_idx,n) = 1;
end

% FH assignment for URLLC UEs
for n=1:Nh
  yi_prime = zeros(Ns,1);

  ue_permuted = randperm(Nu1);
  for i=1:Nu1
    psi_avail = ones(Ns,1) - J(:,n) - yi_prime;
    avail_idx = find(psi_avail>0);
    uei_sc_idx = randperm(length(avail_idx), Ni1);
    psi_ue = zeros(Ns,1);
    psi_ue(avail_idx(uei_sc_idx), 1) = 1;
    yi_prime = yi_prime + psi_ue;

    Y_tmp = Yi{ ue_permuted(i) };
    Y_tmp(:,n) = psi_ue;
    Yi{ ue_permuted(i) } = Y_tmp;
  end
end

JplusI = J;
for i=1:Nu1
  JplusI = JplusI + Yi{i};
end

% FH assignment for mMTC UEs 
%   mMTC UEs are divided into two subgroup in time division

ue_permuted = Nu1 + randperm(Nu2);
for n=1:Nh
  nnn = 1;
  if (mod(n,2))
    nnn = 1;
    ue_permuted = Nu1 + randperm(Nu2);
  else
    nnn = 2;
  end

  yi_prime = zeros(Ns,1);
  %for i=Nu1+nnn : 2 : Nu
  for i=nnn : 2 : Nu2
    psi_avail = ones(Ns,1) - JplusI(:,n) - yi_prime;
    avail_idx = find(psi_avail > 0);
    uei_sc_idx = randperm(length(avail_idx), Ni2);
    psi_ue = zeros(Ns,1);
    psi_ue(avail_idx(uei_sc_idx), 1) = 1;
    yi_prime = yi_prime + psi_ue;

    Y_tmp = Yi{ ue_permuted(i) };
    Y_tmp(:,n) = psi_ue;
    Yi{ ue_permuted(i) } = Y_tmp;
  end
  
end
SNRdB_list = [-2 0 2 4 6 8];
N_realization = 500;
for n=1:N_realization
  fprintf('[D\t n=%d]\n',n);
  for i_snr = 1:length(SNRdB_list)
    SNRdB = SNRdB_list(i_snr);
    output = simQPSK_Y(Yi, J, SNRdB);
    nz_SINR = output.SINR_ues;
    
    y_SINR(i_snr,n) = min(nz_SINR);
    y_SINR_mean(i_snr,n) = mean(nz_SINR);

    output = simBER_givenYiYjam(Yi, J, SNRdB);
    nz_BER = output.BER;
    yDet_BER_min(i_snr,n) = min(nz_BER);
    yDet_BER_p80(i_snr,n) = prctile(nz_BER,80);
    yDet_BER_p90(i_snr,n) = prctile(nz_BER,90);
    yDet_BER_max(i_snr,n) = max(nz_BER);
  end
end
y_SINR = mean(y_SINR,2);
y_SINR_mean = mean(y_SINR_mean,2);
yDet_BER_min = mean(yDet_BER_min,2);
yDet_BER_p80 = mean(yDet_BER_p80,2);
yDet_BER_p90 = mean(yDet_BER_p90,2);
yDet_BER_max = mean(yDet_BER_max,2);

fprintf("output min-SINR for Deterministic Jamming case:\n");
disp(y_SINR)

Yi_det = Yi;
Yjam_det = J;

fprintf("------------- End of Deterministic Case --------------------------\n")
fprintf("------------- ************************* --------------------------\n")
fprintf("------------- Start of Statistical Case --------------------------\n")
fprintf("% Random Jamming FH pattern case\n")
% Random Jamming FH pattern case
PJ = zeros(Ns, Nh);
mid_jb = ceil(Ns/2); 

mu_jam = mid_jb;
sig_jam = Nj_band/ (2^1); % 0.68 interval

pj_obj.a = 1;
pj_obj.b = Ns;
pj_obj.mu = mu_jam;
pj_obj.sigma = sig_jam;
pj_obj.n = Nj;

N_jam_training = 100;
jam_samples=[];
for i=1:N_jam_training
  jam_pattern = genJamPattern(pj_obj);
  jam_samples = [jam_samples; jam_pattern];
end
pj = genPMF_pj( jam_samples, pj_obj.b );
%figure; bar([1:Ns], pj); title('Empirical PMF p_j')

% for f=1:Ns
%   pj(f,1) = 1/(sqrt(2*pi)) / sig_jam * exp(-0.5 * (f-mu_jam)^2 / (sig_jam^2));
% end
for n=1:Nh
    PJ(:,n) = pj;
end


for i=1:Nu
  Yi{i} = zeros(Ns, Nh);
end


% FH assignment for URLLC UEs (Random Jamming Case)
for n=1:Nh
  yi_prime = zeros(Ns,1);
  ue_permuted = randperm(Nu1);
    
  for i=1:Nu1
    p_IJ = (Nj_band * pj + yi_prime)./(Nj_band + ones(1,Ns)*yi_prime);
    psi_i = zeros(Ns,1);

    for ii = 1:Ni1
      [sortVal, sortIdx] = sort(p_IJ); 
      avail_idx = sortIdx( find(sortVal == sortVal(1)) );
      uei_sc_idx = randperm(length(avail_idx), 1);

      psi_i(avail_idx(uei_sc_idx),1) = 1;
      yi_prime_tmp = yi_prime + psi_i;
      p_IJ = (Nj_band * pj + yi_prime_tmp)./(Nj_band + ones(1,Ns)*yi_prime_tmp);
    end

    
    Y_tmp = Yi{ ue_permuted(i) };
    Y_tmp(:,n) = psi_i;
    Yi{ ue_permuted(i) } = Y_tmp;
    yi_prime = yi_prime + psi_i;
  end
end


Interf = zeros(Ns,Nh);
for i=1:Nu1
  Interf = Interf + Yi{i};
end
% FH assignment for mMTC UEs (Random Jammer case)
%   mMTC UEs are divided into two subgroup in time division

ue_permuted = Nu1 + randperm(Nu2);
for n=1:Nh
  nnn = 1;
  if (mod(n,2))
    nnn = 1;
    ue_permuted = Nu1 + randperm(Nu2);
  else
    nnn = 2;
  end

  yi_prime = Interf(:,n);
  for i= nnn : 2 : Nu2
    p_IJ = (Nj_band * pj + yi_prime)./(Nj_band + ones(1,Ns)*yi_prime);
    psi_i = zeros(Ns,1);
    
    for ii = 1:Ni2
      [sortVal, sortIdx] = sort(p_IJ); 
      avail_idx = sortIdx( find(sortVal == sortVal(1)) );
      uei_sc_idx = randperm(length(avail_idx), 1);

      psi_i(avail_idx(uei_sc_idx),1) = 1;
      yi_prime_tmp = yi_prime + psi_i;
      p_IJ = (Nj_band * pj + yi_prime_tmp)./(Nj_band + ones(1,Ns)*yi_prime_tmp);
    end

    
    Y_tmp = Yi{ ue_permuted(i) };
    Y_tmp(:,n) = psi_i;
    Yi{ ue_permuted(i) } = Y_tmp;
    yi_prime = yi_prime + psi_i;
  end
end



% J2 = Nj_band .* PJ;
% generate realization of J indicator matrix
J2 = zeros(Ns, Nh);
for n=1:Nh
  jam_pattern = genJamPattern(pj_obj);
  for ii=1:length(jam_pattern)
    J2(jam_pattern(ii), n) = 1;
  end
end

N_realization = 500;
for n=1:N_realization
  fprintf('[S\t n=%d]\n',n);
  for i_snr = 1:length(SNRdB_list)
    SNRdB = SNRdB_list(i_snr);
    output = simQPSK_Y(Yi, J2, SNRdB);
    nz_SINR = output.SINR_ues;
    
    y2_SINR(i_snr,n) = min(nz_SINR);
    y2_SINR_mean(i_snr,n) = mean(nz_SINR);
  
    output = simBER_givenYiYjam(Yi, J2, SNRdB);
    nz_BER = output.BER;
    yRand_BER_min(i_snr,n) = min(nz_BER);
    yRand_BER_p80(i_snr,n) = prctile(nz_BER,80);
    yRand_BER_p50(i_snr,n) = prctile(nz_BER,50);
    yRand_BER_p30(i_snr,n) = prctile(nz_BER,30);
    yRand_BER_p20(i_snr,n) = prctile(nz_BER,20);
    yRand_BER_p10(i_snr,n) = prctile(nz_BER,10);
    yRand_BER_max(i_snr,n) = max(nz_BER);
  end

end
y2_SINR = mean(y2_SINR,2);
y2_SINR_mean = mean(y2_SINR_mean,2);
yRand_BER_min = mean(yRand_BER_min,2);
yRand_BER_max = mean(yRand_BER_max,2);
yRand_BER_p10 = mean(yRand_BER_p10,2);
yRand_BER_p20 = mean(yRand_BER_p20,2);
yRand_BER_p30 = mean(yRand_BER_p30,2);
yRand_BER_p50 = mean(yRand_BER_p50,2);
yRand_BER_p80 = mean(yRand_BER_p80,2);

fprintf("output min-SINR for Random Jamming case:\n");
disp(y2_SINR)

fprintf("------------- End Data Generation --------------------------\n")
save fh.mat

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

  % for i=1:Nu
  %   Yi_mat = Yi{i};
    
  %   tmpSINR = 0;
  %   count = 0;
  %   for n=1:Nh
  %     yi = Yi_mat(:,n);
  %     y_total = Y_total(:,n);
      
  %     Anz = find(yi > 0);
  %     Dnz = find( (y_total-yi).*yi > 0);
  %     if (length(Anz) > 0)
  %       nhit = length(Dnz) / length(Anz);
  %       tmpSINR = tmpSINR + 10*log10(SNRlin / (1+ SNRlin * nhit));
  %       count = count + 1;
  %     end
  %   end
  %   SINR_ues(i) = tmpSINR / count;
  % end

  output.SINR_ues =SINR_ues; 
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
            % BER(ueI,n) = BER(ueI,n) / (2*Nsym * length(setK));
            BER_array = [BER_array, ber_oneUE/ (2*Nsym*length(setK))];
        end
    end
    BER(:,n) = sort(BER_array);
  end
  output.BER = mean(BER,2);
end