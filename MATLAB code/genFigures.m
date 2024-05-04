close all

set(0, 'defaultlinelinewidth', 1.9);
set(0, 'defaultlinemarkersize', 10);
set(0, 'defaultaxesfontsize', 14);
set(0, 'defaulttextfontsize', 16);
set(0, 'defaultlegendfontsize', 16);

%load fh.mat
% load fh_highSNR.mat
load fh_low2SNR.mat
% load fh_midSNR.mat
% load fh_testSNR

%===========================================================================
% plot SINR vs SNR 
%===========================================================================
figure

plot(SNRdB_list, yNoJam_SINR_p05, '-' );
hold on
plot(SNRdB_list, yDet_noncausal_SINR_p05, 'o' );
plot(SNRdB_list, yStat_m2_SINR_p05, '-s' );
plot(SNRdB_list, yStat_SINR_p05, '^' );
plot(SNRdB_list, yDet_causal_SINR_p05, '--*' );
plot(SNRdB_list, yUniform_SINR_p05, '--' );
% plo'NCKFH',(SNRdB_list, yNoFH_SINR_p05, '-.' );

xlabel('SNR(dB)')
ylabel('SINR(dB)')
legend(...
  'JF', ...
  'NCKFH', ...
  'SFHwP',...
  'SFHwoP',...
  'CKFH',...
  'UFH', ...
  'Location','northwest')
grid on
xlim([-2,7]);  
% ylim([-6,8]);
%mytitleText = [' {\sigma} = ',num2str(sigGau) ];
%title(mytitleText,'Interpreter','tex' );

fig_filename = sprintf("FH_SINR_Ns%d_Nu%d.png", Ns, Nu);
saveas(gcf,fig_filename);

%===========================================================================
% plot BER vs SNR 
%===========================================================================



figure 
% yProbHit = (Nj/Ns) * ones(size(SNRdB_list));
% semilogy(SNRdB_list, yProbHit, 'k--')
semilogy(SNRdB_list([1:5, 8, 10]), yNoJam_BER_p40([1:5, 8, 10]), '-'); 
hold on
semilogy(SNRdB_list([1:5, 8, 10]), yDet_noncausal_BER_p40([1:5, 8, 10]), 'o'); 
semilogy(SNRdB_list([1:5, 8, 10]), yStat_m2_BER_p40([1:5, 8, 10]), '-s'); 
semilogy(SNRdB_list([1:5, 8, 10]), yStat_BER_p40([1:5, 8, 10]), '^'); 
semilogy(SNRdB_list([1:5, 8, 10]), yDet_causal_BER_p40([1:5, 8, 10]), '--*'); 
semilogy(SNRdB_list([1:5, 8, 10]), yUniform_BER_p40([1:5, 8, 10]), '--'); 



SNRlin = 10.^(0.1*SNRdB_list);
BER_theo = 1*qfunc(sqrt(2* SNRlin)); % Theoretical BER 
% semilogy(SNRdB_list, BER_theo)
xlabel('SNR(dB)')
ylabel('BER')
legend(... %     'Prob of Hit', ...
    'JF', ...
    'NCKFH',  ...
    'SFHwP',...
    'SFHwoP',...
    'CKFP',...
    'UFH',  ...
    'Location', 'southwest');
ylim([1e-4,1])
xlim([-1,7]);  
% title(sprintf('Num jammer: %d, Num subcarrier: %d', Nj, Ns))
grid on
fig_filename = sprintf("FH_BER_Ns%d_Nu%d_2.png", Ns, Nu);
saveas(gcf, fig_filename);
fig_filename = sprintf("FH_BER_Ns%d_Nu%d_2.fig", Ns, Nu);
saveas(gcf, fig_filename);

fprintf('size:')
disp(size(nz_BER))
fprintf("------------- Test --------------------------\n")

x=1 - [0.5 1 5 10 20 30 40 50 60 70 80 90 95 99 99.5]/100;
for ii=1:length(SNRdB_list)
    yStat=[yStat_BER_p995(ii), yStat_BER_p99(ii), yStat_BER_p95(ii), yStat_BER_p90(ii), yStat_BER_p80(ii), yStat_BER_p70(ii), ...
        yStat_BER_p60(ii), yStat_BER_p50(ii), yStat_BER_p40(ii), yStat_BER_p30(ii), yStat_BER_p20(ii), yStat_BER_p10(ii),...
        yStat_BER_p05(ii), yStat_BER_p01(ii), yStat_BER_p005(ii) ];
    yStat_m2=[yStat_m2_BER_p995(ii), yStat_m2_BER_p99(ii), yStat_m2_BER_p95(ii), yStat_m2_BER_p90(ii), yStat_m2_BER_p80(ii), yStat_m2_BER_p70(ii), ...
        yStat_m2_BER_p60(ii), yStat_m2_BER_p50(ii), yStat_m2_BER_p40(ii), yStat_m2_BER_p30(ii), yStat_m2_BER_p20(ii), yStat_m2_BER_p10(ii),...
        yStat_m2_BER_p05(ii), yStat_m2_BER_p01(ii), yStat_m2_BER_p005(ii) ];
    yDet_noncausal=[yDet_noncausal_BER_p995(ii), yDet_noncausal_BER_p99(ii), yDet_noncausal_BER_p95(ii), yDet_noncausal_BER_p90(ii), ...
        yDet_noncausal_BER_p80(ii), yDet_noncausal_BER_p70(ii), yDet_noncausal_BER_p60(ii), yDet_noncausal_BER_p50(ii),...
        yDet_noncausal_BER_p40(ii), yDet_noncausal_BER_p30(ii), yDet_noncausal_BER_p20(ii), yDet_noncausal_BER_p10(ii), ...
        yDet_noncausal_BER_p05(ii), yDet_noncausal_BER_p01(ii), yDet_noncausal_BER_p005(ii) ];
    yDet_causal=[yDet_causal_BER_p995(ii), yDet_causal_BER_p99(ii), yDet_causal_BER_p95(ii), yDet_causal_BER_p90(ii), ...
        yDet_causal_BER_p80(ii), yDet_causal_BER_p70(ii), yDet_causal_BER_p60(ii), yDet_causal_BER_p50(ii),...
        yDet_causal_BER_p40(ii), yDet_causal_BER_p30(ii), yDet_causal_BER_p20(ii), yDet_causal_BER_p10(ii), ...
        yDet_causal_BER_p05(ii), yDet_causal_BER_p01(ii), yDet_causal_BER_p005(ii) ];
    yUniform=[yUniform_BER_p995(ii), yUniform_BER_p99(ii), yUniform_BER_p95(ii), yUniform_BER_p90(ii), ...
        yUniform_BER_p80(ii), yUniform_BER_p70(ii), yUniform_BER_p60(ii), yUniform_BER_p50(ii),...
        yUniform_BER_p40(ii), yUniform_BER_p30(ii), yUniform_BER_p20(ii), yUniform_BER_p10(ii), ...
        yUniform_BER_p05(ii), yUniform_BER_p01(ii), yUniform_BER_p005(ii) ];
    yNoJam=[yNoJam_BER_p995(ii), yNoJam_BER_p99(ii), yNoJam_BER_p95(ii), yNoJam_BER_p90(ii), ...
        yNoJam_BER_p80(ii), yNoJam_BER_p70(ii), yNoJam_BER_p60(ii), yNoJam_BER_p50(ii),...
        yNoJam_BER_p40(ii), yNoJam_BER_p30(ii), yNoJam_BER_p20(ii), yNoJam_BER_p10(ii), ...
        yNoJam_BER_p05(ii), yNoJam_BER_p01(ii), yNoJam_BER_p005(ii) ];

    yNoFH=[yNoFH_BER_p995(ii), yNoFH_BER_p99(ii), yNoFH_BER_p95(ii), yNoFH_BER_p90(ii), ...
        yNoFH_BER_p80(ii), yNoFH_BER_p70(ii), yNoFH_BER_p60(ii), yNoFH_BER_p50(ii),...
        yNoFH_BER_p40(ii), yNoFH_BER_p30(ii), yNoFH_BER_p20(ii), yNoFH_BER_p10(ii),...
        yNoFH_BER_p05(ii), yNoFH_BER_p01(ii), yNoFH_BER_p005(ii)];

   
    yProbHit = (Nj/Ns) * ones(size(x));

    figure;
    semilogx(yNoJam,x, '-')
    hold on;
    semilogx(yDet_noncausal,x,'o')
    semilogx(yStat_m2,x,'-s')
    semilogx(yStat,x,'^')
    semilogx(yDet_causal,x,'--*')
    semilogx(yUniform,x,'--')
%     semilogx(yNoFH,x,'-.')
%     semilogx(yProbHit,x, 'k--')

%     title(sprintf('SNR: %d dB', SNRdB_list(ii)))
    ylabel('CDF')
    xlabel('BER')
    legend( 'JF',...
        'NCKFH', ...
        'SFHwP', 'SFHwoP',...
        'CKFH', ...
        'UFH',...
        ... %'NFH',...
        'Location', 'northwest');
    grid on
    xlim([1e-5,1])
    fig_filename = sprintf("FH_BER_percentiel_SNR%d.png", SNRdB_list(ii));
    saveas(gcf, fig_filename);
    fig_filename = sprintf("FH_BER_percentiel_SNR%d.fig", SNRdB_list(ii));
    saveas(gcf, fig_filename);
end



% %===========================================================================
% % plot FH patterns 
% %===========================================================================
% %==============================================
% % (1) plot the FH Determistic Case 
% % (1a) Show interference
% %
% figure
% 
% Yi_plot = 0.5*Yi_det{1};
% for ii=2:Nu
%   Yi_plot = Yi_plot + 0.5 * Yi_det{ii};  
% end
% map = [0,0,0;...
%     0, 0.9, 0.5;...
%     1, 0.5, 0];
% Yi_plot(100,5)=1;
% N=3;
% cmap = colormap(map);
% imagesc(Yi_plot)
% % N=5; cmap = jet(N); colormap(cmap)
% hold on
% L = line(ones(N),ones(N), 'LineWidth',2);
% set(L, {'color'}, mat2cell(cmap, ones(1,N),3));
% legend('empty','UE', 'interference')
% 
% title( sprintf('FH Pattern (I) [Deterministic case]') )
% xlabel('FH Interval')
% ylabel('Subcarrier')
% grid on 
% 
% fig_filename = sprintf("FH_Pattern_det_I.png" );
% saveas(gcf,fig_filename);
% 
% % (1b) Show randomness
% %
% figure
% Yi_plot = Yi_det{1};
% for ii = 2:Nu
%   Yi_plot = Yi_plot - (Yi_plot .* Yi_det{ii}) + ii.* Yi_det{ii};
% end
% imagesc(Yi_plot)
% % 
% % imagesc(Yi_plot, [0,Nu])
% % N=Nu; cmap = jet(N); colormap(cmap)
% 
% title( sprintf('FH Pattern (R) [Deterministic case]') )
% xlabel('FH Interval')
% ylabel('Subcarrier')
% grid on 
% 
% fig_filename = sprintf("FH_Pattern_det_R.png" );
% saveas(gcf,fig_filename);
% 
% 
% 
% 
% figure 
% Yi_plot = Yjam_det;
% % colormap('hot')
% % imagesc(Yi_plot, [0,2])
% % c=colorbar; 
% % c.Ticks = [0 1 2];
% imagesc(Yi_plot)
% N=2; cmap = jet(N); colormap(cmap)
% hold on
% L = line(ones(N),ones(N), 'LineWidth',2);
% set(L, {'color'}, mat2cell(cmap, ones(1,N),3));
% legend('empty','jammer')
% 
% title( sprintf('Jamming Pattern [Deterministic case]') )
% xlabel('FH Interval')
% ylabel('Subcarrier')
% grid on
% 
% fig_filename = sprintf("FH_Jamming_Pattern_det.png" );
% saveas(gcf,fig_filename);
% 
% %============== plot pattern for randomized jamming case ==================
% figure
% Yi_plot = 2*Yi{1};
% for ii=2:Nu
%   Yi_plot = Yi_plot +  (2).* Yi{ii};
%   %Yi_plot = Yi_plot + (ii).* Yi{ii};
%   %fprintf('[U%d, M%d] \n', ii, max(max(Yi_plot(:))) );
% end
% % Yi_plot = Yi_plot + Yjam_det * 257;
% 
% % colormap('hot');
% % imagesc(Yi_plot, [0,2])
% imagesc(Yi_plot)
% N=5; cmap = jet(N); colormap(cmap)
% hold on
% L = line(ones(N),ones(N), 'LineWidth',2);
% set(L, {'color'}, mat2cell(cmap, ones(1,N),3));
% legend('empty','','UE','', 'interference')
% % c=colorbar; 
% % c.Ticks = [0 1 2];
% title( sprintf('FH Pattern (I) [Statistics case]') )
% xlabel('FH Interval')
% ylabel('Subcarrier')
% grid on 
% %colorbar;
% 
% fig_filename = sprintf("FH_Pattern_rand.png" );
% saveas(gcf,fig_filename);
% 
% 
% % (2b) Show randomness
% %
% figure
% Yi_plot = Yi{1};
% for ii = 2:Nu
%   Yi_plot = Yi_plot - (Yi_plot .* Yi_det{ii}) + ii.* Yi_det{ii};
% end
% imagesc(Yi_plot)
% % 
% % imagesc(Yi_plot, [0,Nu])
% % N=Nu; cmap = jet(N); colormap(cmap)
% 
% title( sprintf('FH Pattern (R) [Statistisc case]') )
% xlabel('FH Interval')
% ylabel('Subcarrier')
% grid on 
% 
% fig_filename = sprintf("FH_Pattern_rand_R.png" );
% saveas(gcf,fig_filename);
% 
% 
% figure 
% Yi_plot = J2;
% imagesc(Yi_plot)
% N=2; cmap = jet(N); colormap(cmap)
% hold on
% L = line(ones(N),ones(N), 'LineWidth',2);
% set(L, {'color'}, mat2cell(cmap, ones(1,N),3));
% legend('empty','jammer')
% 
% %c=colorbar; 
% %c.Ticks = [0 1 2];
% title( sprintf('Jamming Pattern [Statistics case]') )
% xlabel('FH Interval')
% ylabel('Subcarrier')
% grid on
% fig_filename = sprintf("FH_Jamming_Pattern_rand.png" );
% saveas(gcf,fig_filename);
% 
% 
% %----------------------
% nz_idx = find(J2>0);
% nz_idx_row = mod( nz_idx-1, Ns) + 1;
% nz_idx_col = ceil( nz_idx / Ns);
% J2_points = [nz_idx_col, nz_idx_row];
% 
% nz_idx = find(Yi{33} > 0);
% nz_idx_row = mod( nz_idx-1, Ns) + 1;
% nz_idx_col = ceil( nz_idx / Ns);
% UE1_points = [nz_idx_col, nz_idx_row];
% 
% nz_idx = find(Yi{34} > 0);
% nz_idx_row = mod( nz_idx-1, Ns) + 1;
% nz_idx_col = ceil( nz_idx / Ns);
% UE2_points = [nz_idx_col, nz_idx_row];
% 
% 
% 
% figure 
% 
% plot(UE1_points(:,1), UE1_points(:,2), 'go', 'MarkerSize', 12)
% hold on
% plot(UE2_points(:,1), UE2_points(:,2), 'bo', 'MarkerSize', 12)
% plot(J2_points(:,1), J2_points(:,2), 'ro', 'MarkerSize', 12)
% 
% 
% % 
% % 
% % set(0, 'defaultlinelinewidth', 5);
% % set(0, 'defaultlinemarkersize', 20);
% figure
% Yi_plot = zeros(Ns,Nh);
% Nuu = 5;
% for ii=1:Nuu
%   Yi_plot = Yi_plot + (ii)/(Nuu+1) * Yi{ii};  
% end
% Yi_plot = Yi_plot + J2;
% 
% map = [1.0,1.0,1.0;...
%     200/255, 200/255, 0.0;...
%     102/255, 204/255, 0.3;...
%     153/255, 51/255, 255/255;...
%     67/255, 181/255, 128/255;...
%     243/255, 213/255, 14/255;...
%     1, 0.1, 0];
% 
% cmap = colormap(map);
% imagesc(Yi_plot)
% % N=5; cmap = jet(N); colormap(cmap)
% N=Nuu+2;
% hold on
% L = line(ones(N),ones(N), 'LineWidth',2);
% set(L, {'color'}, mat2cell(cmap, ones(1,N),3));
% legend(' ','UE1', 'UE2','UE3','UE4','UE5',  ...
%     'Jammer')
% 
% title( sprintf('FH Patterns') )
% xlabel('FH Interval')
% ylabel('Subcarrier')
% grid on 
% 
% fig_filename = sprintf("FH_Pattern_rand_black.png" );
% saveas(gcf,fig_filename);
% 
% fig_filename = sprintf("FH_Pattern_rand_black.fig" );
% saveas(gcf,fig_filename);
% 
