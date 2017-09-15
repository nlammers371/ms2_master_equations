
%%
A = [3,3];
B = [0.3,0.3];
C = [0.03,0.1];
D = [0.05,1];
%%
%Set write paths
subfolder = ['2SpotAuto'];
outpath = ['../../fig/2Spot/' subfolder];
if exist([outpath]) ~= 7
    mkdir(outpath);
end
%%

for a = 1:size(A,2)
Which = A(a);
r1 = B(a);
Kon1 = C(a);
Koff1 = D(a);

vpol = 2; %kb/min
L = 3;%kb
TotalTime = 100; %min
Nreps = 25;
TimeRes = 10; % "
max_lag = 50;
   
ToSave = ['alpha',num2str(Which),'_r',num2str(r1),'_Kon',num2str(Kon1),'_Koff',num2str(Koff1),...
    '_TT',num2str(TotalTime),'_N',num2str(Nreps),'_TRes',num2str(TimeRes),'_vPol',num2str(vpol),'_L',num2str(L)]

alpha = [1,1,1]
alpha(Which) = 2
r2 = alpha(1)*r1;
Kon2 = Koff1/((Koff1+Kon1)/(Kon1*alpha(2))-1);
Koff2 = (Koff1 + (1-alpha(3))*Kon1)/alpha(3);

r = [r1, r2];
Kon = [Kon1, Kon2];
Koff = [Koff1, Koff2];
Ks = [Kon; Koff];
ExpFF = 1+r.*Koff./(Kon+Koff).^2
ExpFFratio = ExpFF(2)/ExpFF(1)

Traces = zeros (TotalTime*60/TimeRes, Nreps);

% mkdir(['/Users/julia/Downloads/Gillespie/',ToSave,'/'])
for n = 1:Nreps

time = []; time(1) = 0;
mRNA = []; mRNA(1) = 0;
State = []; State(1) = 0;
timeInic = [];
K0 = Koff(1); %starts
s=1;
p=1;
Split = TotalTime*60/length(r);
Periods = [1:Split:TotalTime*60];
while time(s) < TotalTime*60
    s = s+1;
    p = floor(time(s-1)./Split)+1;
    dt = 1/K0*log(1/rand(1));
    time(s) = time(s-1) + dt;
    State(s) = ~State(s-1);
    if State(s)==0 & time(s) < TotalTime*60
        cumtimeON = 0;
        while cumtimeON < dt
            dtON = 1/r(p)*log(1/rand(1));
            cumtimeON = cumtimeON + dtON;
            if cumtimeON < dt
                timeInic = [timeInic, time(s-1)+cumtimeON];
            end
        end
    end
    K0 = Ks(State(s)+1,p);
end
State(s+1) = State(s);
time(s+1) = TotalTime*60;
% stairs (time,State-2); hold on
% plot (timeInic, 1-2,'.r');
% plot(Periods,-3,'*k');  hold on
%
%

elT = L / (vpol/60);
TimeScale = [1:TotalTime*60+elT];
Counts = zeros(size(TimeScale));
for f = 2:round(timeInic(end))
    Start = TimeScale(f-1);
    End = TimeScale(f);
    %Counts(f) = sum(timeInic>=Start & timeInic<End)
    Counts(f:(f+elT)) = Counts(f:(f+elT)) + sum(timeInic>=Start & timeInic<End);
end

%Counts(1,end-elT/TimeRes:end) = [];
% plot(TimeScale, Counts,'--b'); hold on
TimeRes = 10; % "
TimeObserved = [1:TimeRes:TotalTime*60];
CountsObserved = Counts(1:TimeRes:TotalTime*60);
% plot(TimeObserved, CountsObserved,'-b'); hold on
% add noise
sigma = 1/10*Kon(1)/(Kon(1)+Koff(1))*r(1)*elT;
noise = normrnd(0, sigma,[1,size(TimeObserved,2)]);
TraceNoise = CountsObserved + noise;
% plot(TimeObserved, TraceNoise); hold on
% hold off
Traces(:,n) = TraceNoise';
% pause(0.05)
% saveas(gcf,['/Users/julia/Downloads/Gillespie/',ToSave,'/',num2str(n),'.pdf'],'pdf')
end
% Mean = figure;
% plot(Traces); hold on
% plot(mean(Traces,2),'.r')
% saveas(Mean,[outpath,ToSave,'_mean.pdf'],'pdf')
% writetable(table(Traces),[outpath,ToSave,'_traces.txt'],'Delimiter','\t')

%
for P = 1:size(Periods,2)
    CumTraces = cumsum(Traces(Periods(P)/TimeRes:(Periods(P)+ ...
        Split)/TimeRes,:),1)
    CumMean(:,P) = mean(CumTraces,2)
    CumVar(:,P) = sum((CumMean(:,P)-CumTraces(:,P)).^2,2)
    FF(:,P) = CumVar(:,P)./CumMean(:,P)
    AutoCorr = zeros(max_lag+1,size(Traces,2));
    for i = 1:Nreps
        AutoCorr(:,i) = autocorr(Traces(ceil(Periods(P)/TimeRes):floor((Periods(P)+ ...
            Split)/TimeRes),i),max_lag);
    end
    ACMean(:,P) = mean(AutoCorr,2)
end
Stats = figure;
plot(CumMean); hold on
plot(CumVar/1000);
plot(FF);
plot(ACMean);
hold off
saveas(Stats,[outpath,ToSave,'_stats.pdf'],'pdf')
FFRatio = FF(:,2)./FF(:,1);
% writetable(table(CumMean,CumVar,FF,FFRatio,ACMean, 'VariableNames',{'CMean','CVar','FF','FFratio','ACMean'}),['/Users/julia/Downloads/Gillespie/',ToSave,'_traces.txt'],'Delimiter','\t')
Ratio = figure;
plot(FFRatio); hold on
plot([1,size(CumMean,1)],[ExpFFratio,ExpFFratio],'b--');
hold off
saveas(Ratio,[outpath,ToSave,'_FFratio.pdf'],'pdf')

Expected(a) = ExpFFratio
Simulation(a) = FFRatio(end)
close all
break
end
%%
All = figure;
plot(Expected,Simulation,'*r'); xlim([0,2]); ylim([0,2]); hold on
plot([0,2],[0,2], 'k--')
saveas(All,[outpath,'expected_vs_simulated.pdf'],'pdf')
