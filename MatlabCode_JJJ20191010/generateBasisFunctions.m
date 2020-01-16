function BF = generateBasisFunctions(Timedata,StartTimes,EndTimes,Options)
%generateBasisFunctions Generate basis functions for lp-ntPET analysis
%   BF=generateBasisFunctions(Timedata,StartTimes,EndTimes,Options)
%   Generates basis functions using Timedata (nx2, double array with frame
%   start and end times), with StartTimes (1xk, double array) describing
%   the onset times of k challenges, EndTimes (1xk, double array)
%   describing the end times of k challenges, and Options containing
%   fields FunctionName (box, exp, gamma), and exp_tau (if exp, 1xj double
%   array) describing exponential decay constant (min), and gamma_alpha (if
%   gamma, 1xj double array describing decay of gamma-variate function.
%
%   Jarkko Johansson, 2018, UMU.se
%

%% Check inputs
errorStat=true;
if size(Timedata,2)==2 && ~isempty(StartTimes) && ~isempty(EndTimes) && ~isempty(Options)
    if size(StartTimes)==size(EndTimes)
        if isstruct(Options)
            if ~isempty(find(contains(fieldnames(Options),'FunctionName')))
                switch Options.FunctionName
                    case 'box'
                        errorStat=false;
                        Func='box';
                    case 'exp'
                        if ~isempty(find(contains(fieldnames(Options),'exp_tau')))
                            errorStat=false;
                            Func='exp';
                            ExpTau=Options.exp_tau;
                        else
                            ErrText=sprintf('FunctionName=%s, but parameter %s not found',Options.FunctionName,'exp_tau');
                        end
                    case 'gamma'
                        if ~isempty(find(contains(fieldnames(Options),'gamma_alpha')))
                            errorStat=false;
                            Func='gamma';
                            GammaAlpha=Options.gamma_alpha;
                        else
                            ErrText=sprintf('FunctionName=%s, but parameter %s not found',Options.FunctionName,'gamma_alpha');
                        end
                    otherwise
                end
            end
        end
    end
end

if ~errorStat
    %fprintf(1,'Ready to go\n');
else
    fprintf(1,'Error in input, %s abort\n',ErrText);
    BF=[];
    return
end

%% Generate Basis-functions
midtimes=mean(Timedata,2);
switch Func
    case 'box'
        for ChalIdx=1:length(StartTimes)
            BFidx=1;
            StartIdx=find(Timedata(:,1)>=StartTimes(ChalIdx),1);
            EndIdx=find(Timedata(:,2)>=EndTimes(ChalIdx),1);
            for Tidx1=StartIdx:EndIdx-1
                for Tidx2=Tidx1+1:EndIdx
                    BF_temp=zeros(length(Timedata),1);
                    BF_temp(Tidx1:Tidx2)=1;
                    BF(:,BFidx,ChalIdx)=BF_temp;
                    BFidx=BFidx+1;
                end
            end
        end
    case 'exp'
        for ChalIdx=1:length(StartTimes)
            BFidx=1;
            StartIdx=find(Timedata(:,1)>=StartTimes(ChalIdx),1);
            EndIdx=find(Timedata(:,2)>=EndTimes(ChalIdx),1);
            for tau=ExpTau
                for Tidx1=StartIdx:EndIdx-1
                    BF_temp=zeros(length(Timedata),1);
                    BF_temp(Tidx1:end)=exp(-tau*(midtimes(Tidx1:end)-midtimes(Tidx1)));
                    BF(:,BFidx,ChalIdx)=BF_temp;
                    BFidx=BFidx+1;
                end
            end
        end
    case 'gamma'
        for ChalIdx=1:length(StartTimes)
            BFidx=1;
            StartIdx=find(Timedata(:,1)>=StartTimes(ChalIdx),1);
            EndIdx=find(Timedata(:,2)>=EndTimes(ChalIdx),1);
            for alpha=GammaAlpha
                for Tidx=StartIdx:EndIdx-1
                    for T2peakidx=Tidx+1:EndIdx
                        BF_temp=zeros(length(Timedata),1);
                        p = [1 alpha midtimes(T2peakidx) midtimes(Tidx)];

                        BF_temp(Tidx:end)=gamma_variate_madsen(p,midtimes(Tidx:end));
                        BF(:,BFidx,ChalIdx)=BF_temp;
                        BFidx=BFidx+1;
                    end
                end
            end
        end
end
end



function y = gamma_variate_madsen(p,t)
%gamma variate formulation of madsen1992 article

gamma=p(1);
alpha=p(2);
tmax=p(3);
delay=p(4);

% idx = find(t > delay);
% dt  = t(idx) - delay;
dt    = t - delay;
dtmax = tmax - delay;

% y(idx) = gamma.*(dt./dtmax).^alpha.*exp(alpha.*((tmax-t(idx))./dtmax));
y = gamma.*(dt./dtmax).^alpha.*exp(alpha.*((tmax-t)./dtmax));
end