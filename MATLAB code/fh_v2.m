clc; clear all;
rng(12345)

Ns = 6;
Nh = 8;
Nu = 3;
Nj = 1;

X_jam = randi(Ns, 1, Nh); % Given hopping pattern of Jammer
X_k_hit = zeros(Nu,Nh);
X_k_succ = zeros(Nu,Nh);
X_k_sinr = zeros(Nu,Nh);

FT_map_hit = zeros(Ns,Nh);
FT_map_succ = zeros(Ns,Nh);
FT_map_sinr = zeros(Ns,Nh);

%fprintf("total combination: %d\n", Ns^Nu);
O1_values=zeros(1,Ns^Nu);
O2_values=zeros(1,Ns^Nu);


%===========================================================================
% Case 1: minimizing hamming distance  
%===========================================================================
for t = 1:Nh
  hit_min = Nu^2 + Nu;
  min_id = -1;

  % Compute at singel Hop, all possible pattern in this hop time
  for pID = 1: Ns^Nu
    
    % Generate one particular pattern for all UEs
    X_set = genPattern(pID, Ns, Nu);

    % Compute number of hits for this pattern
    hit_tmp = 0;
    for ue1 = 1:Nu
      X1 = X_set(ue1,:);

      n_hit =  length(find( X_set - repmat(X1, Nu,1) == 0)) - length(X1);  
      n_hit = n_hit /2 + length(find(X1 == X_jam(t)));
      hit_tmp = hit_tmp + n_hit; 
    end
    
    if (hit_tmp < hit_min)
      hit_min = hit_tmp;
      min_id = pID;
    end
  end
  
  % found optimal pattern at this hop
  % Generate the Freq-Time Occupancy Grid
  for uID = 1:Nu
    pattern_u = mod(floor( (min_id-1)/ (Ns^(uID-1))), Ns) + 1;
    X_k_hit(uID, t) = pattern_u;
    FT_map_hit( pattern_u, t) = FT_map_hit( pattern_u,t) + 1;
  end
  FT_map_hit( X_jam(t), t) = FT_map_hit( X_jam(t), t) + 1;
end

%===========================================================================
% Case 2: maximizing success patterns
%===========================================================================
for t=1:Nh
  succ_max = -1;
  max_id = -1;

  % Compute at singel Hop, all possible pattern in this hop time
  for pID = 1: Ns^Nu
    X_set = genPattern(pID, Ns, Nu);
    
    % Compute number of hits for this pattern
    succ_tmp = 0;
    for f = 1:Ns
      in1 = length(find(X_jam(t) == f));
      in2 = length(find( X_set == f));
      succ_tmp = succ_tmp + length(find( (1+in1)*in2 == 1 ));
    end

    if (succ_tmp > succ_max)
      succ_max = succ_tmp;
      max_id = pID;
    end
  end

  % found optimal pattern
  for uID = 1:Nu
    pattern_u = mod( floor( (max_id-1)/ (Ns^(uID-1))), Ns ) + 1;
    X_k_succ(uID, t) = pattern_u;
    FT_map_succ( pattern_u, t) = FT_map_succ( pattern_u, t) + 1;
  end
  FT_map_succ( X_jam(t), t) = FT_map_succ( X_jam(t), t) + 1;
end

%===========================================================================
% Case 3: maximizing SINR 
%===========================================================================

SNRdB_list = [-2 0 2 4 6 8];
for i_snr = 1:length(SNRdB_list)
  SNRdB = SNRdB_list(i_snr);
  FT_map_sinr = zeros(Ns,Nh);
  for t=1:Nh
    sinr_max = -1;
    max_id = -1;

    % Compute at single Hop, all possible pattern in this hop time
    for pID = 1:Ns^Nu
      X_set = genPattern(pID, Ns, Nu);
      
      % Compute SINR for this pattern
      thisSINR = getSINR(X_set, X_jam(t), SNRdB);
      
      if (thisSINR > sinr_max)
        sinr_max = thisSINR;
        max_id = pID;
      end
    end
    
    % found optimal pattern
    for uID = 1:Nu
      pattern_u = mod(floor( (max_id-1)/ (Ns^(uID-1))), Ns) + 1;
      X_k_sinr(uID, t) = pattern_u;
      FT_map_sinr( pattern_u, t) = FT_map_sinr( pattern_u, t) + 1;
    end
    FT_map_sinr( X_jam(t), t) = FT_map_sinr( X_jam(t), t) + 1;
  end

  SINR_method1(i_snr) = 10 * log10(getSINR(X_k_hit(:,1), X_jam(1), SNRdB)); 
  SINR_method2(i_snr) = 10 * log10(getSINR(X_k_succ(:,1), X_jam(1), SNRdB)); 
  SINR_method3(i_snr) = 10 * log10(getSINR(X_k_sinr(:,1), X_jam(1), SNRdB)); 


  output = simQPSK_MontaCarlo(X_k_hit, X_jam, Ns, SNRdB);
  BER_hit(i_snr) = output.BER;
  output = simQPSK_MontaCarlo(X_k_succ, X_jam, Ns, SNRdB);
  BER_succ(i_snr) = output.BER;
  output = simQPSK_MontaCarlo(X_k_sinr, X_jam, Ns, SNRdB);
  BER_sinr(i_snr) = output.BER;

  fprintf("[SNR:%d]===== max SINR =====\n", SNRdB)
  %fprintf("UE patterns:\n")
  %disp(X_k_sinr)
  fprintf("Freq-Time pattern:\n")
  disp(FT_map_sinr)
end

set(0, 'defaultlinelinewidth', 2);
set(0, 'defaultlinemarkersize', 16);
set(0, 'defaultaxesfontsize', 16);
set(0, 'defaulttextfontsize', 16);
set(0, 'defaultlegendfontsize', 10);
%===========================================================================
% plot SINR vs SNR 
%===========================================================================
figure
plot(SNRdB_list, SINR_method1, '-' );
hold on
plot(SNRdB_list, SINR_method2, '-o' );
plot(SNRdB_list, SINR_method3, '-*' );
xlabel('SNR(dB)')
ylabel('SINR(dB)')
legend('O_{hit}','O_{suc}', 'O_{SINR}', 'Location','northwest')
grid on
%  xlim([0,1]);  
ylim([-6,8]);
%mytitleText = [' {\sigma} = ',num2str(sigGau) ];
%title(mytitleText,'Interpreter','tex' );

fig_filename = sprintf("FH_SINR_Ns%d_Nu%d.png", Ns, Nu);
saveas(gcf,fig_filename);

%===========================================================================
% plot BER vs SNR 
%===========================================================================
  figure(2);
  semilogy(SNRdB_list, BER_hit,'-')  
  hold on
  semilogy(SNRdB_list, BER_succ,'-o')                                 
  semilogy(SNRdB_list, BER_sinr,'-*')                                 
  xlabel('SNR[dB]')                                    
  ylabel('Bit Error Rate');                                         
  legend('O_{hit}','O_{suc}', 'O_{SINR}', 'Location','southwest')
  title(['Probability of Bit Error for QPSK Modulation']);
  grid on;
  ylim([1e-4,1]);
  hold off;

fig_filename = sprintf("FH_BER_Ns%d_Nu%d.png", Ns, Nu);
saveas(gcf,fig_filename);
%===========================================================================
% plot SER vs SNR 
%===========================================================================

%===========================================================================




fprintf("Jammer pattern:\n")
disp(X_jam)

fprintf("===== min number of hit =====\n")
fprintf("UE patterns:\n")
disp(X_k_hit)

fprintf("Freq-Time pattern:\n")
disp(FT_map_hit)


fprintf("===== max success =====\n")
fprintf("UE patterns:\n")
disp(X_k_succ)
fprintf("Freq-Time pattern:\n")
disp(FT_map_succ)

fprintf("===== max SINR =====\n")
fprintf("UE patterns:\n")
disp(X_k_sinr)
fprintf("Freq-Time pattern:\n")
disp(FT_map_sinr)

% fprintf("===== SINR =====\n")
% fprintf("Method 1: min # Hit \n")
% disp(SINR_method1); 
% fprintf("Method 2: max # success \n")
% disp(SINR_method2); 
% fprintf("Method 3: max SINR \n")
% disp(SINR_method3); 
%===========================================================================
function [output] = dist_H(X1, X2)
  output = length( find( X1-X2 ~= 0 ) );
end

function [output] = getSINR(X_set, X_jam, snrdB)
  thisSINR = 1;

  Nu = size(X_set,1);
  Nh = size(X_set,2);

  X_set = [X_set; X_jam];
  snr = 10^(snrdB/10);

  for i=1:Nu
    % each user
    Xi = X_set(i,:);
    tmp = repmat(Xi, size(X_set,1),1) - X_set;
    tmp = tmp(:);

    n_hit = length(find(tmp == 0)) - size(X_set,2);
    thisSINR = thisSINR * snr / ( 1 + n_hit / Nh * snr);
  end
  output = thisSINR^(1/Nu);
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

function [output] = simQPSK_MontaCarlo(X_set, X_jam, Ns, SNRdB)

  Nu = size(X_set, 1);
  Nh = size(X_set, 2);

  N=1e4; % Number of symbols transmited
  b1 = rand(Nu+1, N) > 0.5;
  b2 = rand(Nu+1, N) > 0.5;

  % QPSK sumbol mapping
  I = (2*b1) - 1;
  Q = (2*b2) - 1;
  S = I + 1j*Q;

  SNRlin = 10.^(SNRdB/10);
  N0 = 1./SNRlin;   % Noise Variance

  noise = sqrt(N0/2) * (randn(Ns,N) + 1j*randn(Ns,N)); % AWGN noise (row: subcarrier, col: symbol)
  sig_Rx = zeros(Ns, N) + noise;
  t = 1;

  for ueI = 1:Nu
    % each user, extract its subcarrier index
    k = X_set(ueI,t);
    % superpost the QPSK signals into this subcarrier
    sig_Rx(k,:) = sig_Rx(k,:) + S(ueI,:);
  end

  % add jammer signals
  k = X_jam(t);
  sig_Rx(k,:) = sig_Rx(k,:) + S(Nu+1,:);

  % For BER calculation
  for ueI = 1:Nu
    k = X_set(ueI, t);
    
    sig_I = real(sig_Rx(k,:)); % I component
    sig_Q = imag(sig_Rx(k,:)); % Q component

    bld_I = sig_I > 0; % I decision 
    bld_Q = sig_Q > 0; % Q decision

    b1_error = (bld_I ~= b1(ueI,:)); % Inphase bit error
    b2_error = (bld_Q ~= b2(ueI,:)); % Quadrature bit error

    Error_bit = sum(b1_error(:)) + sum(b2_error(:)); % Total bit error
    BER(ueI,t) = sum(Error_bit)/(2*N); % Simulated BER
    
    % For SER calculation
    error_symbol = or(b1_error(:), b2_error(:)); % if bit in I or bit in Q either wrong than error
    SER(ueI,t) = sum(error_symbol)/N;

    % //BER_theo = 2*qfunc(sqrt(2*SNRlin)); % Theoretical BER 
    % //SER_theo = 2*qfunc(sqrt(2*SNRlin)) - (qfunc(sqrt(2*SNRlin))).^2; % Theoretical SER

  end
  
  


  output.BER = sum(BER(:,t))/Nu;
  output.SER = sum(SER(:,t))/Nu;
end