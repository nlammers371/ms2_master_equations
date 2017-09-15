addpath('../utilities');
% Script to Generate set(s) of traces mimicking 2 spot experiment
% i.e. traces broken into pairs designated as being within same "nucleus"
% Intended to test methdos more separating Extrinsic and Intrinsic Noise
% Terms, and for analyzing each term for information about systems

%--------------------Define Simulation Parameters-------------------------%

%%%General Params

%Specify # Nuclei, average length of observation (seconds) and Tres (seconds)
n_nuclei = 1000;
t_obs = 60*40;
t_res = 20;
avg_seq_length = floor(t_obs / t_res);

%%%Extrinsic Params

%Speficy degree of variabiltiy in trace start and stop times as a percent of t_obs
%. If 0, start and stop times will be synchronous across all traces
% ss_var = .1;
%Specify Variability in initiation rates and transition rates. 
%For now I vary only average rate from nucleus to nucleus. 
%In future temporal fluctuations about spatial mean could be added
r_var = .1;
k_off_var = .1;
k_on_var = .1;

%%%Intrinsic Params

%Set elongation time and determine implied dwell time (note: memory needs to be
%integer for trace generation. t_res and t_elong should be adjusted accordingly)
t_elong =180;
w = round(t_elong / t_res);
%MS2 transcription time in time steps. 
%NL: I've been effectively keep this at zero for the time being
alpha = .001;
%Set num transcription states
K = 2;
%Set noise term as a fraction of mean expression rate
noise = .1;
%Initiation Rates (in AU/s, rather than PolII units)
if K == 2
    r_initiation = [.001, 60];
    %Initial State PDF
    pi0 = [1,0];
elseif K == 3
    r_initiation = [.001, 60, 120];
    %Initial State PDF
    pi0 = [1,0,0];
end
%Degree of correlation between sister chromatids (peg to 1 for now, 
%i.e. no correlation...). Only relevant for 3 state case
beta = 1;
%Number of active_promoters
% n_promoters = 2;
%Elongation rates (bp/s)
% elongation_rates = [2*1000/60];
%If 1, Pol II taken to initiate instantaneously once blocking distance is
%cleared
% instant_fluo = 1;
%Blocking size of Pol II (bp)
% blocking_sizes = [150];

%%%Other Params
%Calibration AU / mRNA
% fluo_per_mRNA = 250;
%Can set multiple values here if you like
k_on_vec = [.05];
k_off_vec = [.5];
%test names
names = {'eve_test'};

%Set write paths
subfolder = ['2SpotTraces_w' num2str(w) '_K' num2str(K)];
outpath = ['../../out/2Spot/' subfolder];
if exist([outpath]) ~= 7
    mkdir(outpath);
end

%----------------------Run Simulations------------------------------------%
%Structure to compile outputs from each simulation
meta_trace_struct = struct;
for j = 1:length(k_off_vec)   
    sim_trace_struct = struct;
    %time vector
    t_vec = 0:t_res:t_obs;
    for i = 1:n_nuclei
        %Draw initiation rates an
        r_nuc = zeros(1,K);
        r_nuc(2) = normrnd(r_initiation(2),r_var*r_initiation(2),1);
        if K == 3
            r_nuc(3) = 2*r_nuc(2);
        end
        %Draw transition rates
        k_on_nuc = normrnd(k_on_vec(j),k_on_vec(j)*k_on_var,1);
        k_off_nuc = normrnd(k_off_vec(j),k_off_vec(j)*k_off_var,1);
        %Define Rate matrix
        if K == 2
            R = [-k_on_nuc, k_off_nuc; 
                  k_on_nuc, -k_off_nuc];
        elseif K == 3
            R = [-2*k_on_nuc, beta*k_off_nuc, 0; 
                  2*k_on_nuc, -beta*(k_off_nuc+k_on_nuc), 2*k_off_nuc;
                  0, beta*k_on_nuc, -2*k_off_nuc];
        end
        
        %Background noise
        nuc_noise = noise*r_nuc(2)*(k_on_nuc/(k_on_nuc+k_off_nuc))*t_res*w;
        for n = 1:2
            seq_length = floor(normrnd(avg_seq_length,ss_var*avg_seq_length,1));
            if seq_length >= length(t_vec)
                seq_length = length(t_vec)-1;
            end
%             trace = synthetic_rate_gillespie_semi_poisson(seq_length, alpha, ...
%                              K, w, R, t_res, r_nuc, nuc_noise, pi0, ...
%                              blocking_sizes, elongation_rates, fluo_per_mRNA, ...
%                              n_promoters, instant_fluo, names, 0);
            trace = synthetic_rate_gillespie(seq_length, alpha, ...
                                K, w, R, t_res, r_nuc, nuc_noise, pi0);

%             trace_sim = trace.loading_scenarios(1);
%             trace_orig = trace.loading_scenarios(end);
            trace_orig = trace.fluo_MS2;seq_length;                        
            ind = (i-1)*2+n;
            sim_trace_struct(ind).fluo = trace_orig;            
            sim_trace_struct(ind).naive_states = trace.naive_states;
            sim_trace_struct(ind).naive_times = trace.transition_times;
            sim_trace_struct(ind).time = t_vec;
            sim_trace_struct(ind).K = K;
            sim_trace_struct(ind).w = w;
            sim_trace_struct(ind).t_res = t_res;
            sim_trace_struct(ind).nucleus = i;
            sim_trace_struct(ind).er = t_elong;
%             sim_trace_struct(ind).blk = trace_sim.blk_time;
%             sim_trace_struct(ind).fluo_per_mRNA = trace_sim.fluo_per_rna;
            sim_trace_struct(ind).nucleus = i;
            sim_trace_struct(ind).r_nuc = r_nuc;
            sim_trace_struct(ind).r_true = r_initiation;
            sim_trace_struct(ind).r_car = r_initiation*r_var;
            sim_trace_struct(ind).R = R;
            sim_trace_struct(ind).k_on_nuc = k_on_nuc;
            sim_trace_struct(ind).k_on_true = k_on_vec(j);            
            sim_trace_struct(ind).k_on_var = k_on_vec(j)*k_on_var;            
            sim_trace_struct(ind).k_off_nuc = k_off_nuc;
            sim_trace_struct(ind).k_off_true = k_off_vec(j);
            sim_trace_struct(ind).k_off_var = k_off_vec(j)*k_off_var;
            sim_trace_struct(ind).noise = nuc_noise;            
        end
    end
    meta_trace_struct(j).simulations = sim_trace_struct;
end
save([outpath '/trace_struct.mat'],'meta_trace_struct');