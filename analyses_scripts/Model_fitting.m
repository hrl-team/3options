%% Simulations for a three-armed bandit, by Sophie Bavard, March 2018

clear

% DATA MARTIX
% COLUMN 1  = participant's number
% COLUMN 2  = phase (0=training, 1=learning, 2=transfer)
% COLUMN 3  = trial number
% COLUMN 4  = condition (negative = forced choice; 0 = transfer)
% COLUMN 5  = left option
% COLUMN 6  = middle option
% COLUMN 7  = right option
% COLUMN 8  = choice (-1=left, 0=middle, 1=righ
% COLUMN 9  = accuracy
% COLUMN 10 = left outcome
% COLUMN 11 = middle outcome
% COLUMN 12 = right outcome
% COLUMN 13 = response time
% COLUMN 14 = participant's code

% EXPLICIT MATRIX
% COLUMN 1  = participant's number
% COLUMN 2  = phase (3=explicit)
% COLUMN 3  = trial number
% COLUMN 4  = option number
% COLUMN 5  = response
% COLUMN 6  = true average value
% COLUMN 7  = participant's code

%% Initialisation

whichmodel = 1:4;

% 1 Q-LEARNING
% 2 DIVISIVE
% 3 RANGE
% 4 RANGE uti
% 5 DIVISIVE uti
% 6 WEBB 2020
% 7 HYBRID D-R
% 8 HYBRID Q-R
% 9 HABIT

phase_fit = listdlg('PromptString','Select a fit',...
    'SelectionMode','single',...
    'ListString',{'Learning','Transfer','Both'});

whichexpe = listdlg('PromptString','Select a version',...
    'SelectionMode','single',...
    'ListString',{'expe1','expe2','expe3'});

priors = listdlg('PromptString','Priors?',...
    'SelectionMode','single',...
    'ListString',{'None','Training'});

if whichexpe==1
    load('data_expe1');
elseif whichexpe==2
    load('data_expe2');
else
    load('data_expe3');
end

subjects = 1:max(data(:,1));

sub=0;
for nsub = subjects

    sub=sub+1;

    % Data matrix
    data_sub = data(data(:,1)==nsub & data(:,2)==1,:) ;

    % Conditions
    condition{sub}  = data(data(:,1)==nsub & data(:,2)==1,4);

    % Learning options and accuracy
    symLT{sub}      = data(data(:,1)==nsub & data(:,2)==1,5:7);
    choLT{sub}      = data(data(:,1)==nsub & data(:,2)==1,9)+1;

    % Learning choice
    if whichexpe < 3
        for i=1:length(data_sub)
            if ~isnan(data_sub(i,8))
                if ismember(symLT{sub}(i,data_sub(i,8)+2), [1 4 6 9])
                    choice{sub}(i,:) = 3;
                elseif ismember(symLT{sub}(i,data_sub(i,8)+2), [2 7])
                    choice{sub}(i,:) = 2;
                elseif ismember(symLT{sub}(i,data_sub(i,8)+2), [3 5 8 10])
                    choice{sub}(i,:) = 1;
                end
            else
                choice{sub}(i,:) = NaN;
            end
        end
    else
        for i=1:length(data_sub)
            if ~isnan(data_sub(i,8))
                if ismember(symLT{sub}(i,data_sub(i,8)+2), [1 4 7 10])
                    choice{sub}(i,:) = 3;
                elseif ismember(symLT{sub}(i,data_sub(i,8)+2), [2 5 8 11])
                    choice{sub}(i,:) = 2;
                elseif ismember(symLT{sub}(i,data_sub(i,8)+2), [3 6 9 12])
                    choice{sub}(i,:) = 1;
                end
            else
                choice{sub}(i,:) = NaN;
            end
        end
    end

    % Learning outcomes
    outcome{sub}    = data(data(:,1)==nsub & data(:,2)==1,10);
    counter1{sub}   = data(data(:,1)==nsub & data(:,2)==1,11);
    counter2{sub}   = data(data(:,1)==nsub & data(:,2)==1,12);

    % Transfer options and choice
    symTT{sub}  = data(data(:,1)==nsub & data(:,2)==2,5:7);
    choTT{sub}  = data(data(:,1)==nsub & data(:,2)==2,8)+2;

    % Transfer outcomes
    Toutcome{sub}    = data(data(:,1)==nsub & data(:,2)==0,10);
    Tcounter1{sub}   = data(data(:,1)==nsub & data(:,2)==0,11);
    Tcounter2{sub}   = data(data(:,1)==nsub & data(:,2)==0,12);

    % Task version
    if whichexpe < 3
        forced(sub) = unique(data(data(:,1)==nsub,15));
    end

    clear data_sub

end

%% OPTIMIZATION

options = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxIter', 10000,'MaxFunEval',10000);
% The option Display is set to off, which means that the optimization algorithm will run silently, without showing the output of each iteration.
% The option MaxIter is set to 10000, which means that the algorithm will perform a maximum of 10,000 iterations.

sub=0;
for nsub = subjects

    sub=sub+1;
    f=waitbar(sub/numel(subjects));

    for model = whichmodel

        if whichexpe < 3

            if model < 4

                params_ini = [1 .5 .5 1 0 0 0 0];
                params_inf = [0 0 0 0 0 0 0 0];
                params_sup = [Inf 1 1 1 0 0 0 0];

            elseif model < 6

                params_ini = [1 .5 .5 1 0 0 0 0];
                params_inf = [0 0 0 0 0 0 0 0];
                params_sup = [Inf 1 1 Inf 0 0 0 0];

            elseif model == 6

                params_ini = [1 .5 .5 1 1 .5 1 00];
                params_inf = [0 0 0 1 0 0 1 0];
                params_sup = [Inf 1 1 1 Inf 1 Inf 0];

            elseif model == 7 || model == 8

                params_ini = [1 .5 .5 .5 0 0 0 0];
                params_inf = [0 0 0 0 0 0 0 0];
                params_sup = [Inf 1 1 1 0 0 0 0];

            elseif model == 9

                params_ini = [1 .5 .5 .5 0 0 0 .5];
                params_inf = [0 0 0 0 0 0 0 0];
                params_sup = [Inf 1 1 1 0 0 0 1];

            end

            [parameters(sub,:,model),ll(sub,model)]=fmincon(@(x) ...
                function_model_simulations_3options(x,Toutcome{sub},Tcounter1{sub},Tcounter2{sub},condition{sub},choice{sub},symTT{sub},choTT{sub},...
                outcome{sub},counter1{sub},counter2{sub},model,1,whichexpe,forced(sub),phase_fit,priors),...
                params_ini,[],[],[],[],params_inf,params_sup,[], options);

        else

            if model < 3

                params_ini = [1 .5 .5 0 0];
                params_inf = [0 0 0 0 0];
                params_sup = [Inf 1 1 0 0];

            elseif model < 5

                params_ini = [1 .5 .5 1 0];
                params_inf = [0 0 0 0 0];
                params_sup = [Inf 1 1 Inf 0];

            else

                params_ini = [1 .5 .5 1 1];
                params_inf = [0 0 0 0 0];
                params_sup = [Inf 1 1 Inf Inf];

            end

            [parameters(sub,:,model),ll(sub,model)]=fmincon(@(x) ...
                function_model_simulations_3options_expe3(x,Toutcome{sub},Tcounter1{sub},Tcounter2{sub},condition{sub},choice{sub},symTT{sub},choTT{sub},...
                outcome{sub},counter1{sub},counter2{sub},model,1,phase_fit,priors),...
                params_ini,[],[],[],[],params_inf,params_sup,[], options);

        end
    end
end

delete(f);

%% BIC comparison

if whichexpe < 3

    nfpm=[3 3 3 4 4 6 4 4 5];

    for n=whichmodel
        bic(:,n)=-2*-ll(:,n)+nfpm(n)*log(180 + 180 * (phase_fit == 3));
    end

end

%% SAVE DATA

if phase_fit == 1
    save(strcat('Optimization_expe',num2str(whichexpe),'_learning'));
elseif fit == 2
    save(strcat('Optimization_expe',num2str(whichexpe),'_transfer'));
elseif phase_fit == 3
    save(strcat('Optimization_expe',num2str(whichexpe),'_both'));
end


