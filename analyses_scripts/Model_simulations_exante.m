%% Simulations for a three-armed bandit, by Sophie Bavard, March 2021

clear

%% Initialisation

whichexpe = listdlg('SelectionMode','single','ListString',{'expe1','expe2'});

% Set up option values for each expe
if whichexpe == 1
    values = [86 50 14 86 14 50 32 14 50 14];
else
    values = [86 50 14 86 14 86 68 50 86 50];
end

% Number of trials
ntrial=45;

% Number of simulated agents
subjects=1:50;

% Task conditions
c1 = zeros(ntrial,numel(subjects));
c2 = zeros(ntrial,numel(subjects));
c3 = zeros(ntrial,numel(subjects));
c4 = zeros(ntrial,numel(subjects));

% Final values
rating_arg  = zeros(10,numel(subjects));

% Choose model to simulate
model = listdlg('PromptString','Select a model',...
    'SelectionMode','single',...
    'ListString',{'Q-LEARNING','DIVISIVE','RANGE'}) ;


%% Simulation

% Number of iterations
repet=100;

% Create the transfer 180x2 matrix (all combinations)
POST = [unique(nchoosek(repmat(1:10,1,2),2),'rows');unique(nchoosek(repmat(1:10,1,2),2),'rows')]; POST(POST(:,1)==POST(:,2),:)=[];

for n = 1:repet
    
    f=waitbar(n/repet);
    
    for sub = 1:numel(subjects)

        % Sample 'random' parameters
        parameters(sub,:) = [gamrnd(1.2,5) betarnd(1.1,1.1) betarnd(1.1,1.1)];
        
        % Create learning matrix (45 trials per condition)
        cond{sub}=[repmat(1,1,ntrial) repmat(2,1,ntrial) repmat(3,1,ntrial) repmat(4,1,ntrial)]'; 
        
        % Randomize trials
        cond{sub} = cond{sub}(randperm(length(cond{sub})));
        
        % Create transfer matrix
        post_sym{sub}=POST;
        
        % Simulation
        [~, proba_correct{sub}, ~, Qfinal(:,sub)] = function_model_simulations_3options(parameters(sub,:),0,0,0,cond{sub},{},post_sym{sub},{},{},{},{},model,2,whichexpe,0,0,0);

        % Simulate transfer with argmax choice rule
        for j = 1:length(post_sym{sub})

            % softmax
%             simchoTT{sub}(j,:)  = 1/(1+exp((Qfinal(post_sym{sub}(j,1)) - Qfinal(post_sym{sub}(j,2))) * parameters(sub,1)));

            % argmax
            if Qfinal(post_sym{sub}(j,1),sub) == Qfinal(post_sym{sub}(j,2),sub)
                simchoTT{sub}(j,:) = rand>0.5;
            else
                simchoTT{sub}(j,:) = max(Qfinal(post_sym{sub}(j,1),sub), Qfinal(post_sym{sub}(j,2),sub)) == Qfinal(post_sym{sub}(j,2),sub);
            end

            % espilon-greedy
%             if rand<0.1
%                 simchoTT{sub}(j,:) = 1 + (rand>0.5);
%             else
%                 simchoTT{sub}(j,:) = 1 + (max(Qfinal(post_sym{sub}(j,1)), Qfinal(post_sym{sub}(j,2))) == Qfinal(post_sym{sub}(j,2)));
%             end

        end

        choice_TT{1,sub} =[simchoTT{sub}(post_sym{sub}(:,2)==1,:);1-simchoTT{sub}(post_sym{sub}(:,1)==1,:)];
        choice_TT{2,sub} =[simchoTT{sub}(post_sym{sub}(:,2)==2,:);1-simchoTT{sub}(post_sym{sub}(:,1)==2,:)];
        choice_TT{3,sub} =[simchoTT{sub}(post_sym{sub}(:,2)==3,:);1-simchoTT{sub}(post_sym{sub}(:,1)==3,:)];
        choice_TT{4,sub} =[simchoTT{sub}(post_sym{sub}(:,2)==4,:);1-simchoTT{sub}(post_sym{sub}(:,1)==4,:)];
        choice_TT{5,sub} =[simchoTT{sub}(post_sym{sub}(:,2)==5,:);1-simchoTT{sub}(post_sym{sub}(:,1)==5,:)];
        choice_TT{6,sub} =[simchoTT{sub}(post_sym{sub}(:,2)==6,:);1-simchoTT{sub}(post_sym{sub}(:,1)==6,:)];
        choice_TT{7,sub} =[simchoTT{sub}(post_sym{sub}(:,2)==7,:);1-simchoTT{sub}(post_sym{sub}(:,1)==7,:)];
        choice_TT{8,sub} =[simchoTT{sub}(post_sym{sub}(:,2)==8,:);1-simchoTT{sub}(post_sym{sub}(:,1)==8,:)];
        choice_TT{9,sub} =[simchoTT{sub}(post_sym{sub}(:,2)==9,:);1-simchoTT{sub}(post_sym{sub}(:,1)==9,:)];
        choice_TT{10,sub}=[simchoTT{sub}(post_sym{sub}(:,2)==10,:);1-simchoTT{sub}(post_sym{sub}(:,1)==10,:)];

        for m=1:10
            rating_arg(m,sub)=rating_arg(m,sub)+nanmean(choice_TT{m,sub})/repet;
        end
        
    end
    
    % Learning test
    
    hyperCond=vector_to_structure_matrix(cond,1,length(cond{1}));
    Proba = vector_to_structure_matrix(proba_correct,1,length(cond{1}));      
    
    c1 = c1 + (structure_matrix_to_plotmatrix(Proba,1,hyperCond,numel(subjects),1,ntrial,0))/repet;
    c2 = c2 + (structure_matrix_to_plotmatrix(Proba,2,hyperCond,numel(subjects),1,ntrial,0))/repet;
    c3 = c3 + (structure_matrix_to_plotmatrix(Proba,3,hyperCond,numel(subjects),1,ntrial,0))/repet;
    c4 = c4 + (structure_matrix_to_plotmatrix(Proba,4,hyperCond,numel(subjects),1,ntrial,0))/repet;
    
end

delete(f);


%% SAVE DATA

if whichexpe==1
    save(strcat('exante_model',num2str(model),'_expe1'), 'rating_arg');
else
    save(strcat('exante_model',num2str(model),'_expe2'), 'rating_arg');
end


