%% THREE ARMED BANDIT


function [lik, Pc3, PcPT, Qf, Qvalues, ws] = function_model_simulations_3options(params,tr,tc1,tc2,s,A,ss,aa,R,CF1,CF2,model,fit,version,forced,phase_fit,priors)

% Parameters
beta    = params(1); % choice temperature
alphaQf = params(2); % alpha update Q-values
alphaQc = params(3); % alpha update counterfactual Q-values

if model == 4 || model == 5
    
    utility = params(4); % curvature
    
end

if model == 6 % Webb 2020
    
    sigma = params(5); % saturation
    w     = params(6); % weight
    pnorm  = params(7); % beta-norm
    
end

if model == 7 || model == 8

    weight = params(4);

end

if model == 9

    weight = params(4);
    hab    = params(8);

    D   = zeros(4,3) + 0.5 + 49.5;
    H   = zeros(4,3) + 0.5;

end



Q   = zeros(4,3) + 0.5 + 49.5 * (model == 1);
L   = zeros(4,3) + 0.5 + 49.5 * (model == 1);

% Qvalues = zeros(4,3) + 0.5 + 49.5 * (model == 1);

r   = zeros(1,length(s));
c1  = zeros(1,length(s));
c2  = zeros(1,length(s));

Rmin = zeros(4,1);%+50;
Rmax = zeros(4,1);%+50;

lik=0;

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
            
            if phase_fit ~= 2 % fit including learning phase
                
                if mod(s(i),2)==1 % conditions with 3 options

                    if ~isnan(a(i)); lik = lik + beta * Q(s(i),a(i)) - log(sum(exp(beta * Q(s(i),:)))); end
                    
                elseif mod(s(i),2)==0 % conditions with 2 options
                    
                    if ~isnan(a(i)); lik = lik + beta * Q(s(i),a(i)) - log(sum(exp(beta * Q(s(i),[1 3])))); end
                    
                end
                
            end
            
        end
        
    elseif fit == 2 % model simulation
        
        
        if mod(s(i),2)==0 % conditions with 2 options
            
            Pc1(i) = 1 / (1 +  exp((Q(s(i),3)-Q(s(i),1))*beta));  % probabilit� de choisir le symbole 1 (le moins bon)
            Pc2(i) =  0  ;  % on "supprime" ce symbole car condition 2 options
            Pc3(i) = 1 / (1 +  exp((Q(s(i),1)-Q(s(i),3))*beta));  % probabilit� de choisir le symbole 3 (le meilleur)
            
        elseif mod(s(i),2)==1 % conditions with 3 options
            
            Pc1(i) = 1 / (1 + exp((Q(s(i),2)-Q(s(i),1))*beta) + exp((Q(s(i),3)-Q(s(i),1))*beta));  % probabilit� de choisir le symbole 1 (le moins bon)
            Pc2(i) = 1 / (1 + exp((Q(s(i),3)-Q(s(i),2))*beta) + exp((Q(s(i),1)-Q(s(i),2))*beta));  % probabilit� de choisir le symbole 2 (le moyen)
            Pc3(i) = 1 / (1 + exp((Q(s(i),1)-Q(s(i),3))*beta) + exp((Q(s(i),2)-Q(s(i),3))*beta));  % probabilit� de choisir le symbole 3 (le meilleur)
            
        end
        
        % make the choice depending on the probabilities to choose each symbol
        
        choice = [1 2 3];                   % symbols
        proba  = [Pc1(i) Pc2(i) Pc3(i)];    % probas
        
        a(i) = choice(find(rand<cumsum(proba),1,'first'));  % choisit un symbole avec sa probabilite
        
    end
    
    % learning steps from here
    
    if ~isnan(a(i))
        
        if s(i) < 0 && forced == 3, complete = 0; else, complete = 1; end
        if s(i) < 0, s(i) = abs(s(i)); end
        
        if fit == 1
            
            r(i) = R(i);
            
            if mod(s(i),2)==0 || a(i)==2
                c1(i) = CF1(i);
                c2(i) = CF2(i);
            else
                c1(i) = CF2(i);
                c2(i) = CF1(i);
            end
            
        else
            
            % for each trial, set the reward behind each symbol
            
            if version == 1
                if s(i)==1 || s(i)==2 % total magnitude
                    outcomes = round(normrnd([14 50 86],4));
                elseif s(i)==3 || s(i)==4 % small magnitude
                    outcomes = round(normrnd([14 32 50],4));
                end
            else
                if s(i)==1 || s(i)==2 % total magnitude
                    outcomes = round(normrnd([14 50 86],4));
                elseif s(i)==3 || s(i)==4 % big magnitude
                    outcomes = round(normrnd([50 68 86],4));
                end
            end
            
            outcomes(outcomes>100) = 100;
            outcomes(outcomes < 0) = 0;
            
            r(i) = outcomes(a(i));
            
            if mod(s(i),2)==0 % conditions with 2 options
                c1(i) = outcomes(4-a(i));
            elseif mod(s(i),2)==1 % conditions with 3 options
                c1(i) = outcomes(mod(a(i)  ,3)+1);
                c2(i) = outcomes(mod(a(i)+1,3)+1);
            end
        end
        
        % update values with reward
        
        if model==1 % Qlearning
            
            deltaR = r(i) - Q(s(i),a(i)) ;
            Q(s(i),a(i)) =  Q(s(i),a(i))   + alphaQf * deltaR;   % chosen option   1 / 2 / 3
            
            if complete
                
                if mod(s(i),2)==0 % conditions with 2 options
                    
                    deltaC = c1(i) - Q(s(i),4-a(i));
                    Q(s(i),4-a(i)) = Q(s(i),4-a(i)) + alphaQc * deltaC;
                    
                elseif mod(s(i),2)==1 % conditions with 3 options
                    
                    deltaC1 = c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
                    deltaC2 = c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
                    Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQc * deltaC1;  % unchosen option 2 / 3 / 1
                    Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQc * deltaC2;  % unchosen option 3 / 1 / 2
                    
                end
            end
            
        elseif model==2 % DIVISIVE normalization (no hidden values)
            
            
            if mod(s(i),2)==0 % conditions with 2 options
                
                if complete
                    opt_sum(i) = (r(i) + c1(i));
                else % when partial feedback, normalize with last seen value
                    opt_sum(i) = (r(i) + L(s(i),4-a(i)));
                end
                
                norm_r(i)  = r(i)  / opt_sum(i);
                norm_c1(i) = c1(i) / opt_sum(i);
                
                deltaR = norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i))   = Q(s(i),a(i))   + alphaQf * deltaR;   % chosen option   1 / 3
                
                if complete
                    deltaC = norm_c1(i) - Q(s(i),4-a(i));
                    Q(s(i),4-a(i)) = Q(s(i),4-a(i)) + alphaQc * deltaC;   % unchosen option  3 / 1
                end
                
            elseif mod(s(i),2)==1 % conditions with 3 options
                
                if complete
                    opt_sum(i) = (r(i) + c1(i) + c2(i));
                else % when partial feedback, normalize with last seen value
                    opt_sum(i) = (r(i) + L(s(i),mod(a(i),3)+1) + L(s(i),mod(a(i)+1,3)+1));
                end
                
                norm_r(i)  =  r(i) / opt_sum(i);
                norm_c1(i) = c1(i) / opt_sum(i);
                norm_c2(i) = c2(i) / opt_sum(i);
                
                deltaR  =  norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i)) = Q(s(i),a(i)) + alphaQf * deltaR;   % chosen option   1 / 2 / 3
                
                if complete
                    deltaC1 = norm_c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
                    deltaC2 = norm_c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
                    Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQc * deltaC1;  % unchosen option 2 / 3 / 1
                    Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQc * deltaC2;  % unchosen option 3 / 1 / 2
                end
            end
            
            L(s(i),    a(i)       ) = r(i);
            L(s(i),mod(a(i)  ,3)+1) = c1(i);
            L(s(i),mod(a(i)+1,3)+1) = c2(i);
            
            
        elseif model==3 % RANGE adaptation (no hidden values)
            
            if mod(s(i),2)==0 % conditions with 2 options
                
                if complete
                    maximum(i) = max([r(i),c1(i)]);
                    minimum(i) = min([r(i),c1(i)]);
                else % when partial feedback, normalize with last seen value
                    maximum(i) = max([r(i),L(s(i),4-a(i))]);
                    minimum(i) = min([r(i),L(s(i),4-a(i))]);
                end
                
                norm_r(i)   =  (r(i)  - minimum(i)) / (maximum(i) - minimum(i));
                norm_c1(i)  =  (c1(i) - minimum(i)) / (maximum(i) - minimum(i));
                
                deltaR = norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i))   = Q(s(i),a(i))   + alphaQf * deltaR;   % chosen option   1 / 2 / 3
                
                if complete
                    deltaC = norm_c1(i) - Q(s(i),4-a(i));
                    Q(s(i),4-a(i)) = Q(s(i),4-a(i)) + alphaQc * deltaC;
                end
                
            elseif mod(s(i),2)==1 % conditions with 3 options
                
                if complete
                    maximum(i) = max([r(i),c1(i),c2(i)]);
                    minimum(i) = min([r(i),c1(i),c2(i)]);
                else % when partial feedback, normalize with last seen value
                    maximum(i) = max([r(i),L(s(i),mod(a(i),3)+1),L(s(i),mod(a(i)+1,3)+1)]);
                    minimum(i) = min([r(i),L(s(i),mod(a(i),3)+1),L(s(i),mod(a(i)+1,3)+1)]);
                end
                
                norm_r(i)   =  (r(i)  - minimum(i)) / (maximum(i) - minimum(i));
                norm_c1(i)  =  (c1(i) - minimum(i)) / (maximum(i) - minimum(i));
                norm_c2(i)  =  (c2(i) - minimum(i)) / (maximum(i) - minimum(i));
                
                deltaR  =  norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i)) = Q(s(i),a(i)) + alphaQf * deltaR;   % chosen option   1 / 2 / 3
                
                if complete
                    deltaC1 = norm_c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
                    deltaC2 = norm_c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
                    Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQc * deltaC1;  % unchosen option 2 / 3 / 1
                    Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQc * deltaC2;  % unchosen option 3 / 1 / 2
                end
            end
            
            L(s(i),    a(i)       ) = r(i);
            L(s(i),mod(a(i)  ,3)+1) = c1(i);
            L(s(i),mod(a(i)+1,3)+1) = c2(i);           
            
            
        elseif model == 4 % RANGE with utility
            
            if mod(s(i),2)==0 % conditions with 2 options
                
                if complete
                    maximum(i) = max([r(i),c1(i)]);
                    minimum(i) = min([r(i),c1(i)]);
                else % when partial feedback, normalize with last seen value
                    maximum(i) = max([r(i),L(s(i),4-a(i))]);
                    minimum(i) = min([r(i),L(s(i),4-a(i))]);
                end
                
                norm_r(i)   =  ((r(i)  - minimum(i)) / (maximum(i) - minimum(i)))^utility;
                norm_c1(i)  =  ((c1(i) - minimum(i)) / (maximum(i) - minimum(i)))^utility;
                
                deltaR = norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i))   = Q(s(i),a(i))   + alphaQf * deltaR;   % chosen option   1 / 2 / 3
                
                if complete
                    deltaC = norm_c1(i) - Q(s(i),4-a(i));
                    Q(s(i),4-a(i)) = Q(s(i),4-a(i)) + alphaQc * deltaC;
                end
                
            elseif mod(s(i),2)==1 % conditions with 3 options
                
                if complete
                    maximum(i) = max([r(i),c1(i),c2(i)]);
                    minimum(i) = min([r(i),c1(i),c2(i)]);
                else % when partial feedback, normalize with last seen value
                    maximum(i) = max([r(i),L(s(i),mod(a(i),3)+1),L(s(i),mod(a(i)+1,3)+1)]);
                    minimum(i) = min([r(i),L(s(i),mod(a(i),3)+1),L(s(i),mod(a(i)+1,3)+1)]);
                end
                
                norm_r(i)   =  ((r(i)  - minimum(i)) / (maximum(i) - minimum(i)))^utility;
                norm_c1(i)  =  ((c1(i) - minimum(i)) / (maximum(i) - minimum(i)))^utility;
                norm_c2(i)  =  ((c2(i) - minimum(i)) / (maximum(i) - minimum(i)))^utility;
                
                deltaR  =  norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i)) = Q(s(i),a(i)) + alphaQf * deltaR;   % chosen option   1 / 2 / 3
                
                if complete
                    deltaC1 = norm_c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
                    deltaC2 = norm_c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
                    Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQc * deltaC1;  % unchosen option 2 / 3 / 1
                    Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQc * deltaC2;  % unchosen option 3 / 1 / 2
                end
            end
            
            L(s(i),    a(i)       ) = r(i);
            L(s(i),mod(a(i)  ,3)+1) = c1(i);
            L(s(i),mod(a(i)+1,3)+1) = c2(i);
            

        elseif model==5 % DIVISIVE with utility
            
            
            if mod(s(i),2)==0 % conditions with 2 options
                
                if complete
                    opt_sum(i) = (r(i) + c1(i));
                else % when partial feedback, normalize with last seen value
                    opt_sum(i) = (r(i) + L(s(i),4-a(i)));
                end
                
                norm_r(i)  = (r(i)  / opt_sum(i))^utility;
                norm_c1(i) = (c1(i) / opt_sum(i))^utility;
                
                deltaR = norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i))   = Q(s(i),a(i))   + alphaQf * deltaR;   % chosen option   1 / 3
                
                if complete
                    deltaC = norm_c1(i) - Q(s(i),4-a(i));
                    Q(s(i),4-a(i)) = Q(s(i),4-a(i)) + alphaQc * deltaC;   % unchosen option  3 / 1
                end
                
            elseif mod(s(i),2)==1 % conditions with 3 options
                
                if complete
                    opt_sum(i) = (r(i) + c1(i) + c2(i));
                else % when partial feedback, normalize with last seen value
                    opt_sum(i) = (r(i) + L(s(i),mod(a(i),3)+1) + L(s(i),mod(a(i)+1,3)+1));
                end
                
                norm_r(i)  = (r(i)  / opt_sum(i))^utility;
                norm_c1(i) = (c1(i) / opt_sum(i))^utility;
                norm_c2(i) = (c2(i) / opt_sum(i))^utility;
                
                deltaR  =  norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i)) = Q(s(i),a(i)) + alphaQf * deltaR;   % chosen option   1 / 2 / 3
                
                if complete
                    deltaC1 = norm_c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
                    deltaC2 = norm_c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
                    Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQc * deltaC1;  % unchosen option 2 / 3 / 1
                    Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQc * deltaC2;  % unchosen option 3 / 1 / 2
                end
            end
            
            L(s(i),    a(i)       ) = r(i);
            L(s(i),mod(a(i)  ,3)+1) = c1(i);
            L(s(i),mod(a(i)+1,3)+1) = c2(i);

                
        elseif model==6 % DIVISIVE normalization (Webb 2020)
            
            
            if mod(s(i),2)==0 % conditions with 2 options
                
                if complete
                    opt_sum(i) = sigma + w * (r(i)^pnorm + c1(i)^pnorm)^(1/pnorm);
                else % when partial feedback, normalize with last seen value
                    opt_sum(i) = sigma + w * (r(i)^pnorm + L(s(i),4-a(i))^pnorm)^(1/pnorm);
                end
                
                norm_r(i)  = r(i)  / opt_sum(i);
                norm_c1(i) = c1(i) / opt_sum(i);
                
                deltaR = norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i))   = Q(s(i),a(i))   + alphaQf * deltaR;   % chosen option   1 / 3
                
                if complete
                    deltaC = norm_c1(i) - Q(s(i),4-a(i));
                    Q(s(i),4-a(i)) = Q(s(i),4-a(i)) + alphaQc * deltaC;   % unchosen option  3 / 1
                end
                
            elseif mod(s(i),2)==1 % conditions with 3 options
                
                if complete
                    opt_sum(i) = sigma + w * (r(i)^pnorm + c1(i)^pnorm + c2(i)^pnorm)^(1/pnorm);
                else % when partial feedback, normalize with last seen value
                    opt_sum(i) = sigma + w * (r(i)^pnorm + L(s(i),mod(a(i),3)+1)^pnorm + L(s(i),mod(a(i)+1,3)+1)^pnorm)^(1/pnorm);
                end
                
                norm_r(i)  =  r(i) / opt_sum(i);
                norm_c1(i) = c1(i) / opt_sum(i);
                norm_c2(i) = c2(i) / opt_sum(i);
                
                deltaR  =  norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i)) = Q(s(i),a(i)) + alphaQf * deltaR;   % chosen option   1 / 2 / 3
                
                if complete
                    deltaC1 = norm_c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
                    deltaC2 = norm_c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
                    Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQc * deltaC1;  % unchosen option 2 / 3 / 1
                    Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQc * deltaC2;  % unchosen option 3 / 1 / 2
                end
            end
            
            L(s(i),    a(i)       ) = r(i);
            L(s(i),mod(a(i)  ,3)+1) = c1(i);
            L(s(i),mod(a(i)+1,3)+1) = c2(i);

        elseif model == 7 % hybrid DIVISIVE-RANGE
            
            if mod(s(i),2)==0 % conditions with 2 options
                
                if complete
                    opt_sum(i) = (r(i) + c1(i));
                    maximum(i) = max([r(i),c1(i)]);
                    minimum(i) = min([r(i),c1(i)]);
                else % when partial feedback, normalize with last seen value
                    opt_sum(i) = (r(i) + L(s(i),4-a(i)));
                    maximum(i) = max([r(i),L(s(i),4-a(i))]);
                    minimum(i) = min([r(i),L(s(i),4-a(i))]);
                end
                
                % DIVISIVE norm
                Dnorm_r(i)  = r(i)  / opt_sum(i);
                Dnorm_c1(i) = c1(i) / opt_sum(i);

                % RANGE norm
                Rnorm_r(i)   =  (r(i)  - minimum(i)) / (maximum(i) - minimum(i));
                Rnorm_c1(i)  =  (c1(i) - minimum(i)) / (maximum(i) - minimum(i));

                % weighted sum
                norm_r(i)  = weight * Rnorm_r(i)  + (1-weight) * Dnorm_r(i);
                norm_c1(i) = weight * Rnorm_c1(i) + (1-weight) * Dnorm_c1(i);

                deltaR = norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i))   = Q(s(i),a(i))   + alphaQf * deltaR;   % chosen option   1 / 3
                
                if complete
                    deltaC = norm_c1(i) - Q(s(i),4-a(i));
                    Q(s(i),4-a(i)) = Q(s(i),4-a(i)) + alphaQc * deltaC;   % unchosen option  3 / 1
                end
                
            elseif mod(s(i),2)==1 % conditions with 3 options
                
                if complete
                    opt_sum(i) = (r(i) + c1(i) + c2(i));
                    maximum(i) = max([r(i),c1(i),c2(i)]);
                    minimum(i) = min([r(i),c1(i),c2(i)]);
                else % when partial feedback, normalize with last seen value
                    opt_sum(i) = (r(i) + L(s(i),mod(a(i),3)+1) + L(s(i),mod(a(i)+1,3)+1));
                    maximum(i) = max([r(i),L(s(i),mod(a(i),3)+1),L(s(i),mod(a(i)+1,3)+1)]);
                    minimum(i) = min([r(i),L(s(i),mod(a(i),3)+1),L(s(i),mod(a(i)+1,3)+1)]);
                end
                
                % DIVISIVE norm
                Dnorm_r(i)  =  r(i) / opt_sum(i);
                Dnorm_c1(i) = c1(i) / opt_sum(i);
                Dnorm_c2(i) = c2(i) / opt_sum(i);
                
                % RANGE norm
                Rnorm_r(i)   =  (r(i)  - minimum(i)) / (maximum(i) - minimum(i));
                Rnorm_c1(i)  =  (c1(i) - minimum(i)) / (maximum(i) - minimum(i));
                Rnorm_c2(i)  =  (c2(i) - minimum(i)) / (maximum(i) - minimum(i));

                % weighted sum
                norm_r(i)  = weight * Rnorm_r(i)  + (1-weight) * Dnorm_r(i);
                norm_c1(i) = weight * Rnorm_c1(i) + (1-weight) * Dnorm_c1(i);
                norm_c2(i) = weight * Rnorm_c2(i) + (1-weight) * Dnorm_c2(i);

                deltaR  =  norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i)) = Q(s(i),a(i)) + alphaQf * deltaR;   % chosen option   1 / 2 / 3
                
                if complete
                    deltaC1 = norm_c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
                    deltaC2 = norm_c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
                    Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQc * deltaC1;  % unchosen option 2 / 3 / 1
                    Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQc * deltaC2;  % unchosen option 3 / 1 / 2
                end
            end

            L(s(i),    a(i)       ) = r(i);
            L(s(i),mod(a(i)  ,3)+1) = c1(i);
            L(s(i),mod(a(i)+1,3)+1) = c2(i);

        elseif model == 8 % hybrid UNBIASED-RANGE


            if mod(s(i),2)==0 % conditions with 2 options

                if complete
                    maximum(i) = max([r(i),c1(i)]);
                    minimum(i) = min([r(i),c1(i)]);
                else % when partial feedback, normalize with last seen value
                    maximum(i) = max([r(i),L(s(i),4-a(i))]);
                    minimum(i) = min([r(i),L(s(i),4-a(i))]);
                end

                % RANGE norm
                Rnorm_r(i)   =  (r(i)  - minimum(i)) / (maximum(i) - minimum(i));
                Rnorm_c1(i)  =  (c1(i) - minimum(i)) / (maximum(i) - minimum(i));

                % weighted sum
                norm_r(i)  = weight * Rnorm_r(i)  + (1-weight) * r(i);
                norm_c1(i) = weight * Rnorm_c1(i) + (1-weight) * c1(i);

                deltaR = norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i))   = Q(s(i),a(i))   + alphaQf * deltaR;   % chosen option   1 / 3
                
                if complete
                    deltaC = norm_c1(i) - Q(s(i),4-a(i));
                    Q(s(i),4-a(i)) = Q(s(i),4-a(i)) + alphaQc * deltaC;   % unchosen option  3 / 1
                end
                
            elseif mod(s(i),2)==1 % conditions with 3 options
                
                if complete
                    maximum(i) = max([r(i),c1(i),c2(i)]);
                    minimum(i) = min([r(i),c1(i),c2(i)]);
                else % when partial feedback, normalize with last seen value
                    maximum(i) = max([r(i),L(s(i),mod(a(i),3)+1),L(s(i),mod(a(i)+1,3)+1)]);
                    minimum(i) = min([r(i),L(s(i),mod(a(i),3)+1),L(s(i),mod(a(i)+1,3)+1)]);
                end
                
                % RANGE norm
                Rnorm_r(i)   =  (r(i)  - minimum(i)) / (maximum(i) - minimum(i));
                Rnorm_c1(i)  =  (c1(i) - minimum(i)) / (maximum(i) - minimum(i));
                Rnorm_c2(i)  =  (c2(i) - minimum(i)) / (maximum(i) - minimum(i));

                % weighted sum
                norm_r(i)  = weight * Rnorm_r(i)  + (1-weight) * r(i);
                norm_c1(i) = weight * Rnorm_c1(i) + (1-weight) * c1(i);
                norm_c2(i) = weight * Rnorm_c2(i) + (1-weight) * c2(i);

                deltaR  =  norm_r(i) - Q(s(i),a(i)) ;
                Q(s(i),a(i)) = Q(s(i),a(i)) + alphaQf * deltaR;   % chosen option   1 / 2 / 3
                
                if complete
                    deltaC1 = norm_c1(i) - Q(s(i),mod(a(i)  ,3)+1) ;
                    deltaC2 = norm_c2(i) - Q(s(i),mod(a(i)+1,3)+1) ;
                    Q(s(i),mod(a(i)  ,3)+1) = Q(s(i),mod(a(i)  ,3)+1)  + alphaQc * deltaC1;  % unchosen option 2 / 3 / 1
                    Q(s(i),mod(a(i)+1,3)+1) = Q(s(i),mod(a(i)+1,3)+1)  + alphaQc * deltaC2;  % unchosen option 3 / 1 / 2
                end
            end
            
            L(s(i),    a(i)       ) = r(i);
            L(s(i),mod(a(i)  ,3)+1) = c1(i);
            L(s(i),mod(a(i)+1,3)+1) = c2(i);

        elseif model == 9 % HABIT

            deltaR = r(i) - D(s(i),a(i)) ;
            D(s(i),a(i)) =  D(s(i),a(i))   + alphaQf * deltaR;   % chosen option   1 / 2 / 3

            if complete

                if mod(s(i),2)==0 % conditions with 2 options

                    deltaC = c1(i) - D(s(i),4-a(i));
                    D(s(i),4-a(i)) = D(s(i),4-a(i)) + alphaQc * deltaC;

                elseif mod(s(i),2)==1 % conditions with 3 options

                    deltaC1 = c1(i) - D(s(i),mod(a(i)  ,3)+1) ;
                    deltaC2 = c2(i) - D(s(i),mod(a(i)+1,3)+1) ;
                    D(s(i),mod(a(i)  ,3)+1) = D(s(i),mod(a(i)  ,3)+1)  + alphaQc * deltaC1;  % unchosen option 2 / 3 / 1
                    D(s(i),mod(a(i)+1,3)+1) = D(s(i),mod(a(i)+1,3)+1)  + alphaQc * deltaC2;  % unchosen option 3 / 1 / 2

                end
            end

            % Habitual controller
            if a(i) == 1
                a_t = [1 0 0];
            elseif a(i) == 2
                a_t = [0 1 0];
            elseif a(i) == 3
                a_t = [0 0 1];
            end

            H(s(i),:)   = H(s(i),:) + hab * (a_t - H(s(i),:));

            % Arbiter
            Q(s(i),:)   = weight * H(s(i),:)  +  (1-weight) * D(s(i),:) ;

        end
    end

    Qvalues(:,:,i) = Q;

end

%% Transfer Test

Qf(1)  = Q(1,3); % 86
Qf(2)  = Q(1,2); % 50
Qf(3)  = Q(1,1); % 14

Qf(4)  = Q(2,3); % 86
Qf(5)  = Q(2,1); % 14

Qf(6)  = Q(3,3); % 86
Qf(7)  = Q(3,2); % 68
Qf(8)  = Q(3,1); % 50

Qf(9)  = Q(4,3); % 86
Qf(10) = Q(4,1); % 50


PcPT=zeros(1,length(ss));

for j = 1:length(ss)
    
    if fit == 1

        if ~isnan(aa(j))
            if (aa(j))

                if phase_fit > 1 % fit including the transfer

                    lik = lik + beta * Qf(ss(j,aa(j))) - log(nansum(exp(beta * Qf(ss(j,~isnan(ss(j,:)))))));

                end

            end
        end
    else
        
        PcPT(j)  = 1/(1+exp((Qf(ss(j,1)) - Qf(ss(j,2))) * beta));

        % make the choice for tranfer        
        choicePT = [1 2];                   % symbols
        probaPT  = [1-PcPT(j) PcPT(j)];     % probas
        
        choPT(j) = choicePT(find(rand<cumsum(probaPT),1,'first'));  % choisit un symbole avec sa probabilite
        
    end
end

if fit == 1
    
    lik = -lik;

    % Prior penalization

    beta    = params(1); % choice temperature
    alphaQf = params(2); % alpha update Q-values
    alphaQc = params(3); % alpha update counterfactual Q-values
    utility = params(4); % curvature
    sigma   = params(5); % saturation
    w       = params(6); % weight
    pnorm   = params(7); % beta-norm

    weight  = params(4); % hybrid model

    hab = params(8); % habit model


    % the parameters are distrubution with mean + variance (different shapes)
    pbeta       = log(gampdf(beta,1.2,5));
    palphaQf    = log(betapdf(alphaQf,1.1,1.1));
    palphaQc    = log(betapdf(alphaQc,1.1,1.1));
    putility    = log(gampdf(utility,1.2,5));
    psigma      = log(gampdf(sigma,2,1));
    pw          = log(gampdf(w,2.1,0.16));
    ppnorm      = log(gampdf(pnorm,2,4));

    pweight     = log(betapdf(weight,1.1,1.1));

    phab    = log(betapdf(hab,1.1,1.1));

    if model == 1       % ABSOLUTE
        p = [pbeta palphaQf palphaQc];
    elseif model == 2   % DIVISIVE
        p = [pbeta palphaQf palphaQc];
    elseif model == 3   % RANGE
        p = [pbeta palphaQf palphaQc];
    elseif model == 4   % RANGE non linear
        p = [pbeta palphaQf palphaQc putility];
    elseif model == 5   % DIVISIVE non linear
        p = [pbeta palphaQf palphaQc putility];
    elseif model == 6   % WEBB 2020
        p = [pbeta palphaQf palphaQc putility psigma pw ppnorm];
    elseif model == 7   % HYBRID
        p = [pbeta palphaQf palphaQc pweight];
    elseif model == 8   % HYBRID
        p = [pbeta palphaQf palphaQc pweight];
    elseif model == 9   % HABIT
        p = [pbeta palphaQf palphaQc pweight phab];
    end

    p = -sum(p);

    lik = p + lik;

end


if fit == 2
    % save structure with necessary variables (model/parameter recovery)
    ws.s = s;
    ws.a = a';
    ws.r = r';
    ws.c1 = c1';
    ws.c2 = c2';
    ws.choPT = choPT';
end

end