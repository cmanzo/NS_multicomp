
% % Copyright (c) 2019-2020 Carlo Manzo & Tina Košuta 

% If you use this code, please cite: 
% (http://dx.doi.org/10.1039/C9CP05616E)
% (https://arxiv.org/abs/1909.13133)

% Permission is hereby granted, free of charge, to any person 
% obtaining a copy of this software and associated documentation files 
% (the "Software"), to deal in the Software without restriction, 
% including without limitation the rights to use, copy, modify, merge, 
% publish, distribute, sublicense, and/or sell copies of the Software, 
% and to permit persons to whom the Software is furnished to do so, 
% subject to the following conditions:

% The above copyright notice and this permission notice shall be 
% included in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%%

clear all
close all
clc

%% set param for Nested Sampling
warning off
% Number of particles
NSpar.N =  30;% 20
% Depth to go to
NSpar.depth =  40;

NSpar.pdfmodel='lognormal'; 
switch NSpar.pdfmodel
    case 'lognormal'
        NSpar.pdf=@(x,mu,sigma) 0.5*(erf( (mu - log(x-1))./(sqrt(2)*sigma) ) -erf( (mu- log(x))./(sqrt(2)*sigma) ) );
        NSpar.pdfpar = [3.3491 0.8462];    
    %    add extra cases
end

NSpar.delt=1.5;

%% load data

locs=readtable('sample_data.txt');
locs=table2array(locs(:,2));
locs=locs(isnan(locs)==0);

%
%% create x
locx=unique(locs);
ll=1:1:3*max(locx); % to make sure of normalization after numerical convolution
n=histc(locs,locx)'; % histogram of localizations

% create calibration function fo monomers
if length(NSpar.pdfpar)==2
    y(:,1)=NSpar.pdf(ll,NSpar.pdfpar(1),NSpar.pdfpar(2));
    % to use one parameter distributions
%elseif length(NSpar.pdfpar)==1 
%    y(:,1)=NSpar.pdf(ll,NSpar.pdfpar(1));
end

% condit to stop
flagmax=0;
% number of species
num_params=1;

% start adding one extra component at the time
while    flagmax<1   % condition on the maximum of logZnew
    num_params = num_params+1;
    
    %% generate random particles (sum = 1)
    particles=dirichletrnd(NSpar.delt*ones(1,num_params),NSpar.N);
    %% add a convolution term
    y(:,num_params)=ifft(fft(y(:,1)).^[num_params]);
    y1=y(locx,:);
    %%
    % calculate log_prior
    logps=logdirichlet(particles, NSpar.delt); % dirichlet
    % sum function (sum(alpha*f)) for each particle
    fmix=particles*y1';  % or mix=y*particles';
    
    % calculate log_likelihood on every particle
    %logls= sum((n).*log(fmix),2);
    logls=(n*log(fmix'))';
    % initialize NS parameters and evidence
    logwidth=0;
    logZnew=0;
    %%
    keep=[];
    keep2=[];
    % initialize ratio of evidence for convergence
    Zrat=Inf;
    i=0;
    while Zrat>log(1E-5)
        i=i+1;
        % Find worst particle
        [worst_logl ,worst] = min(logls);
        %calculate logweight times loglikelihood
        logwidth=-i*log(1+1/NSpar.N) -log(NSpar.N) +worst_logl  ;
        
        % dead particles
        keep = [keep; particles(worst,:) logls(worst)];
        keep2=[keep2; logwidth];
        
        % set loglikelihood threshold
        threshold = logls(worst);
        
        % Copy a surviving particles
        flag=0;
        while flag==0
            which = randsample(NSpar.N, 1);
            if(which ~= worst)
                break
            end
        end
        particles(worst,:) = particles(which,:);
        % upgrade likelihood and prior 
        logps(worst) = logps(which);
        logls(worst) = logls(which);
        
        % evolve particle with MCMC
        sig=.1;% tune step width based on rejection 
        accepted=0; % initialize accepted
        R=0;
        for j=1:NSpar.depth
            new0=particles(worst,:);
 
             kk=randsample(num_params, 1);
             step_size=random('norm',0,sig,1);
%             % periodic boundaries in [0,1]             
             new0(1,kk)=new0(1,kk) + step_size;
             new0(1,kk)=new0(1,kk)-floor(new0(1,kk));       
             new0(1,:)=new0(1,:)./sum(new0(1,:));

                 
                fmix_new=new0*y1';  % or mix=y*particles';
                logl_new=(n*log(fmix_new'))';% log likelihood
                logp_new=logdirichlet(new0(1,:), NSpar.delt);% prior
                
                % Metropolis-Hastings
                if logl_new > threshold && rand<=exp(min([logp_new-logps(worst),0]))
                    particles(worst,:) = new0(1,:);
                    logls(worst) = logl_new;
                    logps(worst) = logp_new;
                    accepted = accepted + 1;
                else
                  R=R+1;  
                end
                
                if accepted>R
                    sig=sig*exp(1/accepted);
                elseif accepted<R
                    sig=sig/exp(1/R);
                end
              
            %end
        end %end MCMC
        clear fmix_new logl_new logp_new 
        
        logZnew=(log(sum( exp( keep2-max(keep2) )) ) + max(keep2));
        logZres=-i*log(1+1/NSpar.N) -log(NSpar.N) +  (log(sum( exp( logls-max(logls) )) ) + max(logls)); 
        
        
        Zrat=(logZres   -  logZnew);
        
    end
    
    %%
    clear Zrat fmix logZnew logZres logwidth new0 z z0 worst worst_logl 
    
    logwidth2=-i*log(1+1/NSpar.N) -log(NSpar.N)  + logls; %  
    
    % add dead particles
    keep2=[keep2; logwidth2];
    keep=[keep; particles logls];
    
    clear particles logls logwidth2 logps sig 
    
    % calculate Evidence
    NSout.logZnew(num_params-1)=(log(sum( exp( keep2-max(keep2) )) ) + max(keep2));
    
   % posterior weights (normalized) 
    post_weights = exp( keep2   - NSout.logZnew(num_params-1));
    
    
    % entropy
    NSout.ent(num_params-1) = -sum(post_weights.*log(post_weights + 1E-300));
    
    % effective posterior sample size
    NSout.ess(num_params-1) = floor(exp(NSout.ent(num_params-1)));
    
    % Information
    NSout.H(num_params-1) =  sum(post_weights.*(keep(:, size(keep,2)) -  NSout.logZnew(num_params-1)));
    
    % Error of logEvidence
    NSout.err(num_params-1) = sqrt(NSout.H(num_params-1)./NSpar.N);
    if NSout.H(num_params-1)<0
        return
    end
    
    
    %% weighted mean
    NSout.alphas{num_params-1}=sum(keep(:,1:end-1).*post_weights);
    NSout.erralpha{num_params-1}= sqrt(sum( (keep(:,1:end-1) -NSout.alphas{num_params-1}).^2 .*post_weights));
            
    % check for max conditions to stop model
    [m1,in1]=max(NSout.logZnew);

    if num_params-1>2
        if in1<num_params-1
            flagmax=1;
        end
    end
    clear m1 in1
    
end
%% display parameters corresponding to max logZnew
[NSout.logZ,imaxZ]=max(NSout.logZnew);
%
NSout_r.NbestZ=imaxZ+1;
% 
NSout_r.alphas=NSout.alphas{imaxZ};
%
NSout_r.erralpha=NSout.erralpha{imaxZ};
%
NSout_r.err=NSout.err(imaxZ);
%
NSout_r.ent=NSout.ent(imaxZ);
NSout_r.ess=NSout.ess(imaxZ);
NSout_r.H=NSout.H(imaxZ);
NSout_r
 




