addpath('./utilities');
% Script to Test Analytic Results for Intrinsic Noise Given System Params
%--------------------Define Simulation Parameters-------------------------%
%Set write paths

% subfolder = '3state_nn_cont_section5_1';
subfolder = '3state_nn_val';
outpath = ['../figs/fano_figs/' subfolder];
if exist([outpath]) ~= 7
    mkdir(outpath);
end
out_name = 'nn';
%Specify # Data Points and T res to simulate
T = 80*60;
deltaT = 10;
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
    r_emission = [.001, 60];
    %Initial State PDF
    pi0 = [0,1];
elseif K == 3
    r_emission = [.001, 60, 120];
    %Initial State PDF
    pi0 = [0,0,1];
end
%Calibration AU / mRNA
fluo_per_rna = 300;
%Obsolete Variable
wait_time = 10;
%Number of traces to generate
n_traces = 500;
% Param details in title
param_details = 0;
% Plot continuous
plot_continuous = 0;
% Plot nn case?
plot_nn = 0;
% k_on_vec = .005:.05:.55;
% k_off_vec = fliplr(.005:.05:.55);

k_on_vec = [.008];
k_off_vec = fliplr([.014]);
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
    %Noise due to background fluctuations
    noise = .1*w*r_emission(2)*deltaT*k_on_vec(j)/(k_on_vec(j)+k_off_vec(j));
    %-----------------------Run Simulation------------------------------------%
    trace_mat = zeros(seq_length,n_traces);
    cm_trace_mat = zeros(seq_length,n_traces);
    cm_trace_mat_nn = zeros(seq_length,n_traces);
    cm_trace_mat_cont = zeros(seq_length,n_traces);
    %Iterate
    on_steps = [];
    for i = 1:n_traces
        trace = synthetic_rate_gillespie_two(seq_length, alpha, ...
                                        K, w, R, deltaT, r_emission, noise, pi0, ...
                                        fluo_per_rna, wait_time, conv);
        trace_mat(:,i) = trace.fluo_MS2;
        cm_trace_mat(:,i) = cumsum(trace.fluo_MS2);
        cm_trace_mat_nn(:,i) = cumsum(trace.fluo_MS2_no_noise);
        cm_trace_mat_cont(:,i) = cumsum(trace.fluo_MS2_orig_noise);
        on_steps = [on_steps sum(trace.naive_states==2)];
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
    
    var_predicted_p = factor*((2*(avg_rate^2)*k_on*k_off)/(k_on+k_off)^3 + ...
        (k_on*avg_rate)/(k_on+k_off))*time_vec;
    
    var_predicted_c = factor*((2*(avg_rate^2)*k_on*k_off)/(k_on+k_off)^3  ...
        )*time_vec;
    
    var_vec_p = var(cm_trace_mat'/(fluo_per_rna*w))' ;
    var_vec_c = var(cm_trace_mat_cont'/(fluo_per_rna*w))' ;
    var_vec_nn = var(cm_trace_mat_nn'/(fluo_per_rna*w))' ;
    
    plot(var_vec_p, 'Color', cm(1,:),'LineWidth',1.5);
    plot(var_predicted_p,'LineStyle','--','Color', cm(30,:),'LineWidth',1.5);
    
    if plot_continuous
        plot(var_vec_c, 'Color', cm(45,:),'LineWidth',1.5);
        plot(var_predicted_c,'LineStyle','--','Color', cm(60,:),'LineWidth',1.5);
        legend('Poisson', 'Poisson predicted', 'Continuous', ...
            'Continuous predicted','Location','southeast')
        
    elseif plot_nn 
        plot(var_vec_nn, 'Color', cm(45,:),'LineWidth',1.5);
        legend('Simulated (Noise)', 'Predicted', 'Simulated (No Noise)', ...
            'Location','southeast')
    
    else
        legend('Simulated', 'Predicted','Location','southeast')
    end
    if param_details
        title(strvcat('                        Variance: Predicted vs. Actual',...
            ['k_{on}:' num2str(k_on) '| k_{off}:' num2str(k_off), '| r:' num2str(r_emission(2)) 'AU | FluoPerRNA:' num2str(fluo_per_rna) ' | NTraces:' ...
            num2str(n_traces) ' | Noise:' num2str(noise) 'AU | dT:' num2str(deltaT)]), 'fontsize',14);
    else
        vt = title(['Variance: Predicted vs. Actual (' num2str(K) ' States)']);
    end
    xlabel('Time Steps');
    saveas(var_fig, [outpath '/' 'var_fig_' num2str(j) out_name '.png'],'png');
    
    %%% Get Mean %%% 
    mean_fig = figure('Visible','off');
    hold on
    avg_rate = r_emission(2) / fluo_per_rna;
    time_vec = 1:seq_length;
    time_vec = time_vec.*deltaT;
    mean_predicted = factor*avg_rate*(k_on/(k_on + k_off))*time_vec;

    mean_vec_p = mean(cm_trace_mat,2) / w / fluo_per_rna;
    mean_vec_c = mean(cm_trace_mat_cont,2) / w / fluo_per_rna; 
    plot(mean_vec_p,'Color', cm(1,:),'LineWidth',1.5);
    plot(mean_predicted,'LineStyle','--','Color' ,cm(20,:),'LineWidth',1.5);
    if plot_continuous
        plot(mean_vec_c,'Color', cm(45,:),'LineWidth',1.5)
        legend('Poisson','Continuous', 'Predicted (Both)','Location','southeast');
    else
        legend('Simulated','Predicted','Location','southeast');
    end
    if param_details
        title(strvcat('                         Mean: Predicted vs. Actual',...
            ['k_{on}:' num2str(k_on) '| k_{off}:' num2str(k_off), '| r:' num2str(r_emission(2)) 'AU | FluoPerRNA:' num2str(fluo_per_rna) ' | NTraces:' ...
            num2str(n_traces) ' | Noise:' num2str(noise) 'AU | dT:' num2str(deltaT)]), 'fontsize',14);
    else
        mt = title(['Mean: Predicted vs. Actual (' num2str(K) ' States)']);
    end
    xlabel('Time Steps');
    ylabel('AU');
    saveas(mean_fig, [outpath '/' 'mean_fig_' num2str(j) out_name '.png'],'png');
    
    %%% Get  Fano %%%
    fano_fig = figure('Visible','off');
    hold on
    fano_predicted_p = var_predicted_p / mean_predicted;
    fano_predicted_c = var_predicted_c / mean_predicted;
    
    fano_vec_p = var_vec_p ./ mean_vec_p;
    fano_vec_c = var_vec_c ./ mean_vec_c;
    
    plot(fano_vec_p,'Color', cm(1,:),'LineWidth',1.5);
    plot(1:length(fano_vec_p),linspace(fano_predicted_p,fano_predicted_p,length(fano_vec_p)),'--','Color', cm(30,:),'LineWidth',1.5);
    if plot_continuous
        plot(fano_vec_c,'LineStyle','--','Color', cm(45,:),'LineWidth',1.5)
        plot(1:length(fano_vec_c),linspace(fano_predicted_c,fano_predicted_c,length(fano_vec_c)),'--','Color', cm(60,:),'LineWidth',1.5);
        legend('Poisson','Poisson Predicted', 'Continuous', 'Continuous predicted','Location','southeast');
    else
        legend('Simulated','Predicted','Location','southeast');
    end
    if param_details
        title(strvcat('                     Fano Factor: Predicted vs. Actual',...
            ['k_{on}:' num2str(k_on) '| k_{off}:' num2str(k_off), '| r:' num2str(r_emission(2)) 'AU | FluoPerRNA:' num2str(fluo_per_rna) ' | NTraces:' ...
            num2str(n_traces) ' | Noise:' num2str(noise) 'AU | dT:' num2str(deltaT)]), 'fontsize',14);
    else
        title(['Fano Factor: Predicted vs. Actual (' num2str(K) 'States)']);
    end
    xlabel('Time Steps');
    saveas(fano_fig, [outpath '/' 'fano_fig_' num2str(j) out_name '.png'],'png');
    
    trifecta = figure('Visible','off');
    trifecta.PaperPosition = [0 0 18 5];
    subplot(1,3,1);
    hold on
    plot(mean_vec_p, 'Color', cm(1,:),'LineWidth',1.5);
    plot(mean_predicted,'LineStyle','--','Color', cm(30,:),'LineWidth',1.5);
    if plot_continuous
        plot(mean_predicted,'LineStyle','--','Color', cm(45,:),'LineWidth',1.5);
        plot(mean_vec_c,'Color', cm(45,:),'LineWidth',1.5)
    end
    title(['Mean: Predicted vs. Actual (' num2str(K) ' States)'], 'fontsize',14);
    xlabel('Time Steps');
    
    subplot(1,3,2);
    hold on
    plot(var_vec_p, 'Color', cm(1,:),'LineWidth',1.5);
    plot(var_predicted_p,'LineStyle','--','Color', cm(30,:),'LineWidth',1.5);
    if plot_continuous
        plot(var_vec_c, 'Color', cm(45,:),'LineWidth',1.5);
        plot(var_predicted_c,'LineStyle','--','Color', cm(60,:),'LineWidth',1.5);
    end
    title(['Variance: Predicted vs. Actual (' num2str(K) ' States)'],  'fontsize',14);
    xlabel('Time Steps');
    
    subplot(1,3,3);
    hold on
    plot(fano_vec_p,'Color', cm(1,:),'LineWidth',1.5);
    plot(1:length(fano_vec_p),linspace(fano_predicted_p,fano_predicted_p,length(fano_vec_p)),'--','Color', cm(30,:),'LineWidth',1.5);
    if plot_continuous
        plot(fano_vec_c, 'Color', cm(45,:),'LineWidth',1.5);
        plot(1:length(fano_vec_c),linspace(fano_predicted_c,fano_predicted_c,length(fano_vec_c)),'LineStyle','--','Color', cm(60,:),'LineWidth',1.5);
    end
    title(['Fano Factor: Predicted vs. Actual (' num2str(K) ' States)'],  'fontsize',14);
    if plot_continuous
        legend('Simulated Poisson', 'Predicted Poisson','Simulated Continuous', 'Predicted Continuous', 'Location' , 'southeast');
    else
        legend('Simulated', 'Predicted', 'Location' , 'southeast');
    end
    xlabel('Time Steps');
    saveas(trifecta, [outpath '/' 'full_fig_' num2str(j) '_' out_name '.png'],'png');
end

% Make Trace Figures
% t_vec = (1:seq_length) .*deltaT;
% %Raw Traces
% tr_fig = figure;
% hold on
% plot(repmat(t_vec,n_traces,1)', trace_mat ./ fluo_per_rna,'LineWidth',.5);
% 
% title('Simulated Traces');
% xlabel('Time (Time Steps)');
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
% xlabel('Time (Time Steps)');
% ylabel('N mRNA');
% grid on 
% saveas(cm_fig, [outpath '/cm_activity.eps'], 'epsc');
% saveas(cm_fig, [outpath '/cm_activity.png'], 'png');
% 

% delete(pool)    