%% THREE ARMED BANDIT


function [lik, Pc3, PcPT, Qf] = function_model_simulations_3options_expe3(params,tr,tc1,tc2,s,A,ss,aa,R,CF1,CF2,model,fit,phase_fit,priors)

% Parameters
beta    = params(1); % choice temperature
alphaQf = params(2); % alpha update Q-values
alphaQc = params(3); % alpha update counterfactual Q-values
omegaQf = params(4);
omegaQc = params(5);

Q   = zeros(4,3) + 0.5 + 49.5 * (model == 1);
r   = zeros(1,length(s));
c1  = zeros(1,length(s));
c2  = zeros(1,length(s));

lik=0;
x = 0:1:100;

if priors == 2
    for t = 1:length(tr)
        Rmin = Rmin + min([tr(t) tc1(t) tc2(t)])/length(tr);
        Rmax = Rmax + max([tr(t) tc1(t) tc2(t)])/length(tr);
    end
end

for i = 1:length(s)
    
    % softmax: update probabilities to choose each symbol with Qvalues
    
    if fit == 1 % model fitting
        
        a(i) = A(i);
        
        if s(i)>0 % free choice
            
            if ismember(phase_fit,[1 3]) % fit including learning phase
                
                if ~isnan(a(i)); lik = lik + beta * Q(s(i),a(i)) - log(sum(exp(beta * Q(s(i),:)))); end
                
            end
            
        else % forced choice between options 1 and 2
            
            if a(i)==3, disp('error'); end
            
            if ismember(phase_fit,[1 3]) % fit including learning phase
                
                if ~isnan(a(i)); lik = lik + beta * Q(abs(s(i)),a(i)) - log(sum(exp(beta * Q(abs(s(i)),1:2)))); end
                
            end
            
        end
        
    elseif fit == 2 % model simulation
        
        if s(i)>0 % conditions with free choices only
            
            Pc1(i) = 1 / (1 + exp((Q(s(i),2)-Q(s(i),1))*beta) + exp((Q(s(i),3)-Q(s(i),1))*beta));  % probability de choisir le symbole 1 (le moins bon)
            Pc2(i) = 1 / (1 + exp((Q(s(i),3)-Q(s(i),2))*beta) + exp((Q(s(i),1)-Q(s(i),2))*beta));  % probability de choisir le symbole 2 (le moyen)
            Pc3(i) = 1 / (1 + exp((Q(s(i),1)-Q(s(i),3))*beta) + exp((Q(s(i),2)-Q(s(i),3))*beta));  % probability de choisir le symbole 3 (le meilleur)
            
        elseif s(i)<0 % conditions with forced choices
            
            Pc1(i) =  1 / (1 +  exp((Q(abs(s(i)),2)-Q(abs(s(i)),1))*beta));
            Pc2(i) =  1 / (1 +  exp((Q(abs(s(i)),1)-Q(abs(s(i)),2))*beta));
            Pc3(i) =  0;
            
        end
        
        % make the choice depending on the probabilities to choose each symbol
        
        choice = [1 2 3];                   % symbols
        proba  = [Pc1(i) Pc2(i) Pc3(i)];    % probas
        
        a(i) = choice(find(rand<cumsum(proba),1,'first'));  % choisit un symbole avec sa probabilite
        
    end
    
    % learning steps from here
    
    if ~isnan(a(i))
        
        if s(i) < 0, s(i) = abs(s(i)); end
        
        if fit == 1
            
            r(i) = R(i);
            
            if a(i)==2
                c1(i) = CF1(i);
                c2(i) = CF2(i);
            else
                c1(i) = CF2(i);
                c2(i) = CF1(i);
            end
            
        else
            
            % for each trial, set the reward behind each symbol
            
            if s(i)==1 || s(i)==2 % total magnitude
                outcomes = round(normrnd([14 50 86],4));
            elseif s(i)==3 || s(i)==4 % small magnitude
                outcomes = round(normrnd([14 32 50],4));
            end
            
            outcomes(outcomes>100) = 100;
            outcomes(outcomes < 0) = 0;
            
            r(i) = outcomes(a(i));
            
            c1(i) = outcomes(mod(a(i)  ,3)+1);
            c2(i) = outcomes(mod(a(i)+1,3)+1);

        end
        
        % update values with reward
        
        if model==1 % Qlearning 2 alphas
            
            deltaI = r(i) - Q(s(i),a(i)) ;
            Q(s(i),a(i)) =  Q(s(i),a(i))   + alphaQf * deltaI;   % chosen option   1 / 2 / 3
            
            deltaC1 = c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
            deltaC2 = c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
            Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQc * deltaC1;  % unchosen option 2 / 3 / 1
            Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQc * deltaC2;  % unchosen option 3 / 1 / 2
                       
            
        elseif model==2 % RANGE 2 alphas
            
            norm_r(i)  =  (r(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)]));
            norm_c1(i) = (c1(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)]));
            norm_c2(i) = (c2(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)]));
            
            deltaI  =  norm_r(i) - Q(s(i),a(i)) ;
            deltaC1 = norm_c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
            deltaC2 = norm_c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
            
            Q(s(i),a(i))            = Q(s(i),a(i))             + alphaQf * deltaI;   % chosen option   1 / 2 / 3
            
            Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQc * deltaC1;  % unchosen option 2 / 3 / 1
            Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQc * deltaC2;  % unchosen option 3 / 1 / 2
                      
            
        elseif model==3 % RANGE 1 omega
            
            norm_r(i)  =  ((r(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)])))^omegaQf;
            norm_c1(i) = ((c1(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)])))^omegaQf;
            norm_c2(i) = ((c2(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)])))^omegaQf;
            
            deltaI  =  norm_r(i) - Q(s(i),a(i)) ;
            deltaC1 = norm_c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
            deltaC2 = norm_c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
            
            Q(s(i),a(i))            = Q(s(i),a(i))             + alphaQf * deltaI;   % chosen option   1 / 2 / 3
            
            Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQc * deltaC1;  % unchosen option 2 / 3 / 1
            Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQc * deltaC2;  % unchosen option 3 / 1 / 2
                      
            
        elseif model==4 % RANGE 1 omega 1 alpha
            
            norm_r(i)  =  ((r(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)])))^omegaQf;
            norm_c1(i) = ((c1(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)])))^omegaQf;
            norm_c2(i) = ((c2(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)])))^omegaQf;
            
            deltaI  =  norm_r(i) - Q(s(i),a(i)) ;
            deltaC1 = norm_c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
            deltaC2 = norm_c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
            
            Q(s(i),a(i))            = Q(s(i),a(i))             + alphaQf * deltaI;   % chosen option   1 / 2 / 3
            
            Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQf * deltaC1;  % unchosen option 2 / 3 / 1
            Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQf * deltaC2;  % unchosen option 3 / 1 / 2

        elseif model==5 % RANGE 2 omega
            
            norm_r(i)  =  ((r(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)])))^omegaQf;
            norm_c1(i) = ((c1(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)])))^omegaQc;
            norm_c2(i) = ((c2(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)])))^omegaQc;
            
            deltaI  =  norm_r(i) - Q(s(i),a(i)) ;
            deltaC1 = norm_c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
            deltaC2 = norm_c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
            
            Q(s(i),a(i))            = Q(s(i),a(i))             + alphaQf * deltaI;   % chosen option   1 / 2 / 3
            
            Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQc * deltaC1;  % unchosen option 2 / 3 / 1
            Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQc * deltaC2;  % unchosen option 3 / 1 / 2

        elseif model==6 % RANGE 2 omega 1 alpha
            
            norm_r(i)  =  ((r(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)])))^omegaQf;
            norm_c1(i) = ((c1(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)])))^omegaQc;
            norm_c2(i) = ((c2(i) - min([r(i),c1(i),c2(i)])) / (max([r(i),c1(i),c2(i)]) - min([r(i),c1(i),c2(i)])))^omegaQc;
            
            deltaI  =  norm_r(i) - Q(s(i),a(i)) ;
            deltaC1 = norm_c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
            deltaC2 = norm_c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
            
            Q(s(i),a(i))            = Q(s(i),a(i))             + alphaQf * deltaI;   % chosen option   1 / 2 / 3
            
            Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQf * deltaC1;  % unchosen option 2 / 3 / 1
            Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQf * deltaC2;  % unchosen option 3 / 1 / 2
            
        end
    end
end


%% Transfer Test

Qf(1)  = Q(1,3); % 86
Qf(2)  = Q(1,2); % 50
Qf(3)  = Q(1,1); % 14

Qf(4)  = Q(2,3); % 86
Qf(5)  = Q(2,2); % 50
Qf(6)  = Q(2,1); % 14

Qf(7)  = Q(3,3); % 50
Qf(8)  = Q(3,2); % 32
Qf(9)  = Q(3,1); % 14

Qf(10) = Q(4,3); % 50
Qf(11) = Q(4,2); % 32
Qf(12) = Q(4,1); % 14

PcPT=zeros(1,length(ss));

for j = 1:length(ss)
    
    if fit == 1
        
        if ~isnan(aa(j))
            if (aa(j))
                
                if ismember(phase_fit,[2 3]) % fit including the transfer
                    
                    lik = lik + beta * Qf(ss(j,aa(j))) - log(nansum(exp(beta * Qf(ss(j,~isnan(ss(j,:)))))));
                    
                end
                
            end
        end
    else
        
        PcPT(j)  = 1/(1+exp((Qf(ss(j,1)) - Qf(ss(j,2))) * beta));
        
    end
end

%%

lik=-lik;

% Prior penalization

if fit==1

    beta    = params(1); % choice temperature
    alphaQf = params(2); % alpha update Q-values
    alphaQc = params(3); % alpha update counterfactual Q-values
    omegaQf = params(4);
    omegaQc = params(5);

    % the parameters are distrubution with mean + variance (different shapes)
    pbeta     = log(gampdf(beta,1.2,5));
    palphaQf  = log(betapdf(alphaQf,1.1,1.1));
    palphaQc  = log(betapdf(alphaQc,1.1,1.1));
    pomegaQf  = log(gampdf(omegaQf,1.2,5));
    pomegaQc  = log(gampdf(omegaQc,1.2,5));

    if model == 1
        p = [pbeta palphaQf palphaQc];
    elseif model == 2
        p = [pbeta palphaQf palphaQc];
    elseif model == 3
        p = [pbeta palphaQf palphaQc pomegaQf];
    elseif model == 4
        p = [pbeta palphaQf pomegaQf];
    elseif model == 5
        p = [pbeta palphaQf palphaQc pomegaQf pomegaQc];
    elseif model == 6
        p = [pbeta palphaQf pomegaQf pomegaQc];
    end

    p = -sum(p);

    lik = p + lik;

end

end