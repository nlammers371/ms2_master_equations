addpath('./utilities');
% Script to Test Analytic Results for Intrinsic Noise Given System Params
%--------------------Define Simulation Parameters-------------------------%
%Set write paths

subfolder = '3state_no_corr';
outpath = ['../figs/fano_figs/' subfolder];
if exist([outpath]) ~= 7
    mkdir(outpath);
end
out_name = 'no_corr';
%Specify # Data Points and T res to simulate
T = 40*240;
deltaT = 18.75;
seq_length = floor(T / deltaT);

%Set elongation time and determine implied memory
t_elong = 150;
w = round(t_elong / deltaT);
%MS2 transcription time (t steps)
alpha = 0;
%Set num transcription states
K = 3;
%Correlation Factor (beta>1: postive corr, <1, anti-corr)
beta = 1;
%Matrix to convert compound state to individual promoter states
if K == 2
    conv = [1, 1; 1, 2];
elseif K == 3
    conv = [1, 1; 1, 2; 2,2];
end
%Emission Rates
if K == 2
    r_emission = [.001, 30];
    %Initial State PDF
    pi0 = [0,1];
elseif K == 3
    r_emission = [.001, 30, 60];
    %Initial State PDF
    pi0 = [0,0,1];
end
%Noise due to background fluctuations
noise = 0;
%Calibration AU / mRNA
fluo_per_rna = 50;
%Obsolete Variable
wait_time = 10;
%Number of traces to generate
n_traces = 500;
%Poisson or Continuous Generator?
Poisson = 1;

k_on_vec = .005:.045:.5;
k_off_vec = fliplr(.005:.045:.5);

% k_on_vec = [.008];
% k_off_vec = fliplr([.014]);
% pool = parpool(4);
for j = 1:length(k_off_vec)
    %Define Rate matrix
    if K == 2
        R = [-k_on_vec(j), k_off_vec(j); 
              k_on_vec(j), -k_off_vec(j)];
    elseif K == 3
        R = [-2*k_on_vec(j), beta*k_off_vec(j), 0; 
              2*k_on_vec(j), -beta*(k_off_vec(j)+k_on_vec(j)), 2*k_off_vec(j);
              0, beta*k_on_vec(j), -2*k_off_vec(j)];
    end
    %-----------------------Run Simulation------------------------------------%
    trace_mat = zeros(seq_length,n_traces);
    cm_trace_mat = zeros(seq_length,n_traces);
    cm_trace_mat_cont = zeros(seq_length,n_traces);
    %Iterate
    for i = 1:n_traces
        trace = synthetic_rate_gillespie_two(seq_length, alpha, ...
                                        K, w, R, deltaT, r_emission, noise, pi0, ...
                                        fluo_per_rna, wait_time, conv);
        trace_mat(:,i) = trace.fluo_MS2;
        cm_trace_mat(:,i) = cumsum(trace.fluo_MS2);
        cm_trace_mat_cont(:,i) = cumsum(trace.fluo_MS2_orig_noise);
    end
    colormap('winter');
    cm = colormap;
    % Get  Variance 
    var_fig = figure('Visible','off');
    hold on
    avg_rate = r_emission(2) / fluo_per_rna;
    time_vec = 1:seq_length;
    time_vec = time_vec(1:end).*deltaT;
    if K == 2
        factor = 1;
        k_on = R(2,1);
        k_off = R(1,2);
    elseif K == 3
        factor = 2;
        k_on = R(3,2);
        k_off = R(1,2);
    end
    
    var_prediction = factor*((2*(avg_rate^2)*k_on*k_off)/(k_on+k_off)^3 + ...
        (k_on*avg_rate)/(k_on+k_off))*time_vec;
    
    var_vec = var(cm_trace_mat'/(fluo_per_rna*w))' ;
    var_vec_cont = var(cm_trace_mat_cont'/(fluo_per_rna*w))' ;
    plot(var_vec, 'Color', cm(5,:),'LineWidth',1.5);
    plot(var_vec_cont, 'Color', cm(40,:),'LineWidth',1.5);
    plot(var_prediction,'LineStyle','--','Color', cm(55,:),'LineWidth',1.5);
    legend('Poisson','Continuous', 'Prediction','Location','southeast')
    title(strvcat('                     Variance: Predicted vs. Actual',...
        ['k_{on}:' num2str(k_on) '| k_{off}:' num2str(k_off), '| r:' num2str(r_emission(2)) 'AU | FluoPerRNA:' num2str(fluo_per_rna) ' | NTraces:' ...
        num2str(n_traces) ' | Noise:' num2str(noise) 'AU']), 'fontsize',10);
    saveas(var_fig, [outpath '/' 'var_fig_' num2str(j) out_name '.png'],'png');
    
    %%% Get Mean %%% 
    mean_fig = figure('Visible','off');
    hold on
    avg_rate = r_emission(2) / fluo_per_rna;
    time_vec = 1:seq_length;
    time_vec = time_vec.*deltaT;
    mean_prediction = factor*avg_rate*(k_on/(k_on + k_off))*time_vec;

    mean_vec = mean(cm_trace_mat,2) / w / fluo_per_rna;
    mean_vec_cont = mean(cm_trace_mat_cont,2) / w / fluo_per_rna; 
    plot(mean_vec,'Color', cm(5,:),'LineWidth',1.5)
    plot(mean_vec_cont,'Color', cm(40,:),'LineWidth',1.5)
    plot(mean_prediction,'LineStyle','--','Color' ,cm(55,:),'LineWidth',1.5);
    legend('Poisson','Continuous', 'Prediction','Location','southeast');
     title(strvcat('                      Mean: Predicted vs. Actual',...
        ['k_{on}:' num2str(k_on) '| k_{off}:' num2str(k_off), '| r:' num2str(r_emission(2)) 'AU | FluoPerRNA:' num2str(fluo_per_rna) ' | NTraces:' ...
        num2str(n_traces) ' | Noise:' num2str(noise) 'AU']), 'fontsize',10);
    saveas(mean_fig, [outpath '/' 'mean_fig_' num2str(j) out_name '.png'],'png');
    
    %%% Get  Fano %%%
    fano_fig = figure('Visible','off');
    hold on
    fano_prediction = var_prediction / mean_prediction;

    fano_vec = var_vec ./ mean_vec;
    fano_vec_cont = var_vec_cont ./ mean_vec_cont;
    plot(fano_vec,'Color', cm(5,:),'LineWidth',1.5)
    plot(fano_vec_cont,'LineStyle','--','Color', cm(40,:),'LineWidth',1.5)
    plot(1:length(fano_vec),linspace(fano_prediction,fano_prediction,length(fano_vec)),'--','Color', cm(55,:),'LineWidth',1.5);
    legend('Poisson','Continuous', 'Prediction','Location','southeast');
     title(strvcat('                 Fano Factor: Predicted vs. Actual',...
        ['k_{on}:' num2str(k_on) '| k_{off}:' num2str(k_off), '| r:' num2str(r_emission(2)) 'AU | FluoPerRNA:' num2str(fluo_per_rna) ' | NTraces:' ...
        num2str(n_traces) ' | Noise:' num2str(noise) 'AU']), 'fontsize',10);
    saveas(fano_fig, [outpath '/' 'fano_fig_' num2str(j) out_name '.png'],'png');
end

% Make Trace Figures
% t_vec = (1:seq_length) .*deltaT;
% %Raw Traces
% tr_fig = figure;
% hold on
% plot(repmat(t_vec,n_traces,1)', trace_mat ./ fluo_per_rna,'LineWidth',.5);
% 
% title('Simulated Traces');
% xlabel('Time (seconds)');
% ylabel('N Active PolII');
% grid on 
% % saveas(tr_fig, [outpath '/raw_traces.eps'], 'epsc');
% % saveas(tr_fig, [outpath '/raw_traces.png'], 'png');
% 
% %Raw Traces
% cm_fig = figure;
% hold on
% plot(repmat(t_vec,n_traces,1)', cm_trace_mat ./ fluo_per_rna,'LineWidth',.5);
% 
% title('Cumulative mRNA Produced by Simulated Traces');
% xlabel('Time (seconds)');
% ylabel('N mRNA');
% grid on 
% saveas(cm_fig, [outpath '/cm_activity.eps'], 'epsc');
% saveas(cm_fig, [outpath '/cm_activity.png'], 'png');
% 

% delete(pool)    