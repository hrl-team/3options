%% Simulations for a three-armed bandit, by Sophie Bavard, March 2021

clear

%% Initialisation

phase_fit = listdlg('PromptString','Select a version',...
    'SelectionMode','single',...
    'ListString',{'Learning','Transfer','Both'});

whichexpe = listdlg('PromptString','Select a version',...
    'SelectionMode','single',...
    'ListString',{'expe1','expe2','expe3'});

if phase_fit == 1
    load(strcat('Optimization_expe',num2str(whichexpe),'_learning'));
elseif phase_fit == 2
    load(strcat('Optimization_expe',num2str(whichexpe),'_transfer'));
elseif phase_fit == 3
    load(strcat('Optimization_expe',num2str(whichexpe),'_both'));
end

% Set up option values for each expe
if whichexpe == 1
    values = [86 50 14 86 14 50 32 14 50 14];
elseif whichexpe == 2
    values = [86 50 14 86 14 86 68 50 86 50];
else
    values = [86 50 14 86 50 14 50 32 14 50 32 14];
end

% Number of trials
ntrial=45;

% Task conditions
c1 = zeros(ntrial,numel(subjects));
c2 = zeros(ntrial,numel(subjects));
c3 = zeros(ntrial,numel(subjects));
c4 = zeros(ntrial,numel(subjects));

% Final values
rating_m=zeros((10+2*(whichexpe==3)),numel(subjects));

% Choose model to simulate
if whichexpe < 3
    model = listdlg('PromptString','Select a model',...
        'SelectionMode','single',...
        'ListString',{'Q-LEARNING','DIVISIVE','RANGE','RANGE uti','DIVISIVE uti / RANGE w+ (expe3)','DIVISIVE Webb','HYBRID D-R','HYBRID Q-R','HABIT'}) ;
else
    model = listdlg('PromptString','Select a model',...
        'SelectionMode','single',...
        'ListString',{'Q-LEARNING','RANGE','RANGE 1 omega 2 alpha','RANGE 1 omega 1 alpha','RANGE 2 omega 2 alpha','RANGE 2 omega 1 alpha'}) ;
end

parameters = parameters(:,:,model);

%% Simulation

% Number of iterations
repet=100;

% Create the transfer 180x2 or 132x2 matrix (all combinations)
if whichexpe < 3
    POST = [unique(nchoosek(repmat(1:10,1,2),2),'rows');unique(nchoosek(repmat(1:10,1,2),2),'rows')]; POST(POST(:,1)==POST(:,2),:)=[];
else
    POST = [unique(nchoosek(repmat(1:12,1,1),2),'rows');unique(nchoosek(repmat(1:12,1,1),2),'rows')]; POST(POST(:,1)==POST(:,2),:)=[];
end

for n = 1:repet
    
    f=waitbar(n/repet);
    
    for sub = 1:numel(subjects)
        
        % Create learning matrix (45 trials per condition)
        if whichexpe < 3
            cond{sub}=[repmat(1,1,ntrial) repmat(2,1,ntrial) repmat(3,1,ntrial) repmat(4,1,ntrial)]';
        else
            if sub < 51
                cond{sub}=[repmat(1,1,ntrial) repmat(2,1,25) repmat(-2,1,20) repmat(3,1,ntrial) repmat(4,1,25) repmat(-4,1,20)]'; % 20 forced trials
            elseif sub < 101
                cond{sub}=[repmat(1,1,ntrial) repmat(2,1,20) repmat(-2,1,25) repmat(3,1,ntrial) repmat(4,1,20) repmat(-4,1,25)]'; % 25 forced trials
            elseif sub < 151
                cond{sub}=[repmat(1,1,ntrial) repmat(2,1,12) repmat(-2,1,33) repmat(3,1,ntrial) repmat(4,1,12) repmat(-4,1,33)]'; % 33 forced trials
            else
                cond{sub}=[repmat(1,1,ntrial) repmat(2,1,11) repmat(-2,1,34) repmat(3,1,ntrial) repmat(4,1,11) repmat(-4,1,34)]'; % 34 forced trials
            end
        end

        % Randomize trials
        cond{sub} = cond{sub}(randperm(length(cond{sub})));
        
        % Create transfer matrix
        post_sym{sub}=POST;
        
        % Simulation
        if whichexpe < 3
            [~, proba_correct{sub}, proba_transfer{sub}, Qfinal(:,sub)] = function_model_simulations_3options(parameters(sub,:),{},{},{},cond{sub},{},post_sym{sub},{},{},{},{},model,2,whichexpe,forced(sub),0,0);
        else
            [~, proba_correct{sub}, proba_transfer{sub}, Qfinal(:,sub)] = function_model_simulations_3options_expe3(parameters(sub,:),{},{},{},cond{sub},{},post_sym{sub},{},{},{},{},model,2,0,0);
        end

        proba_transfer{sub} = proba_transfer{sub}';
        cond{sub} = abs(cond{sub});

        % Transfer test choice rate
        
        choice_TT{1,sub} =[proba_transfer{sub}(post_sym{sub}(:,2)==1,:);1-proba_transfer{sub}(post_sym{sub}(:,1)==1,:)];
        choice_TT{2,sub} =[proba_transfer{sub}(post_sym{sub}(:,2)==2,:);1-proba_transfer{sub}(post_sym{sub}(:,1)==2,:)];
        choice_TT{3,sub} =[proba_transfer{sub}(post_sym{sub}(:,2)==3,:);1-proba_transfer{sub}(post_sym{sub}(:,1)==3,:)];
        choice_TT{4,sub} =[proba_transfer{sub}(post_sym{sub}(:,2)==4,:);1-proba_transfer{sub}(post_sym{sub}(:,1)==4,:)];
        choice_TT{5,sub} =[proba_transfer{sub}(post_sym{sub}(:,2)==5,:);1-proba_transfer{sub}(post_sym{sub}(:,1)==5,:)];
        choice_TT{6,sub} =[proba_transfer{sub}(post_sym{sub}(:,2)==6,:);1-proba_transfer{sub}(post_sym{sub}(:,1)==6,:)];
        choice_TT{7,sub} =[proba_transfer{sub}(post_sym{sub}(:,2)==7,:);1-proba_transfer{sub}(post_sym{sub}(:,1)==7,:)];
        choice_TT{8,sub} =[proba_transfer{sub}(post_sym{sub}(:,2)==8,:);1-proba_transfer{sub}(post_sym{sub}(:,1)==8,:)];
        choice_TT{9,sub} =[proba_transfer{sub}(post_sym{sub}(:,2)==9,:);1-proba_transfer{sub}(post_sym{sub}(:,1)==9,:)];
        choice_TT{10,sub}=[proba_transfer{sub}(post_sym{sub}(:,2)==10,:);1-proba_transfer{sub}(post_sym{sub}(:,1)==10,:)];
        if whichexpe==3
            choice_TT{11,sub} =[proba_transfer{sub}(post_sym{sub}(:,2)==11,:);1-proba_transfer{sub}(post_sym{sub}(:,1)==11,:)];
            choice_TT{12,sub} =[proba_transfer{sub}(post_sym{sub}(:,2)==12,:);1-proba_transfer{sub}(post_sym{sub}(:,1)==12,:)];
        end

        for m=1:(10+2*(whichexpe==3))
            rating_m(m,sub)=rating_m(m,sub)+nanmean(choice_TT{m,sub})/repet;
        end
        
        perfPT_m(:,sub) = nanmean([proba_transfer{sub}(values(post_sym{sub}(:,2))>values(post_sym{sub}(:,1)));...
            1-proba_transfer{sub}(values(post_sym{sub}(:,2))<values(post_sym{sub}(:,1)))]);
        
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

perf_m = [nanmean(c1); nanmean(c2); nanmean(c3); nanmean(c4)];

%% Retrieve data

sub=0;
for nsub = subjects
    
    sub=sub+1;

    if whichexpe < 3

        cond1(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)==1,9);
        cond2(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)==2,9);
        cond3(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)==3,9);
        cond4(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)==4,9);

        cho{sub}  = data(data(:,1)==nsub & data(:,2)==2,8);
        sym{sub}  = data(data(:,1)==nsub & data(:,2)==2,5:7);

        rating(1 ,sub) = sum((sym{sub}(:,1)==1 & cho{sub}==-1) + (sym{sub}(:,2)==1 & cho{sub}==0) + (sym{sub}(:,3)==1 & cho{sub}==1))/36;
        rating(2 ,sub) = sum((sym{sub}(:,1)==2 & cho{sub}==-1) + (sym{sub}(:,2)==2 & cho{sub}==0) + (sym{sub}(:,3)==2 & cho{sub}==1))/36;
        rating(3 ,sub) = sum((sym{sub}(:,1)==3 & cho{sub}==-1) + (sym{sub}(:,2)==3 & cho{sub}==0) + (sym{sub}(:,3)==3 & cho{sub}==1))/36;
        rating(4 ,sub) = sum((sym{sub}(:,1)==4 & cho{sub}==-1) + (sym{sub}(:,2)==4 & cho{sub}==0) + (sym{sub}(:,3)==4 & cho{sub}==1))/36;
        rating(5 ,sub) = sum((sym{sub}(:,1)==5 & cho{sub}==-1) + (sym{sub}(:,2)==5 & cho{sub}==0) + (sym{sub}(:,3)==5 & cho{sub}==1))/36;
        rating(6 ,sub) = sum((sym{sub}(:,1)==6 & cho{sub}==-1) + (sym{sub}(:,2)==6 & cho{sub}==0) + (sym{sub}(:,3)==6 & cho{sub}==1))/36;
        rating(7 ,sub) = sum((sym{sub}(:,1)==7 & cho{sub}==-1) + (sym{sub}(:,2)==7 & cho{sub}==0) + (sym{sub}(:,3)==7 & cho{sub}==1))/36;
        rating(8 ,sub) = sum((sym{sub}(:,1)==8 & cho{sub}==-1) + (sym{sub}(:,2)==8 & cho{sub}==0) + (sym{sub}(:,3)==8 & cho{sub}==1))/36;
        rating(9 ,sub) = sum((sym{sub}(:,1)==9 & cho{sub}==-1) + (sym{sub}(:,2)==9 & cho{sub}==0) + (sym{sub}(:,3)==9 & cho{sub}==1))/36;
        rating(10,sub) = sum((sym{sub}(:,1)==10& cho{sub}==-1) + (sym{sub}(:,2)==10& cho{sub}==0) + (sym{sub}(:,3)==10& cho{sub}==1))/36;

    else

        cond1(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)==1,9);
        cond2(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & abs(data(:,4))==2,9);
        cond3(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)==3,9);
        cond4(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & abs(data(:,4))==4,9);

        sym{sub}  = data(data(:,1)==nsub & data(:,2)==2,5:7);
        cho{sub}  = data(data(:,1)==nsub & data(:,2)==2,8);

        for i=1:12
            rating(i,sub) = sum((sym{sub}(:,1)==i & cho{sub}==-1) + (sym{sub}(:,2)==i & cho{sub}==0) + (sym{sub}(:,3)==i & cho{sub}==1))/22;
        end

    end
end

%% SAVE DATA

if phase_fit == 1
    save(strcat('simulation_model',num2str(model),'_expe',num2str(whichexpe),'_learning'));
elseif phase_fit == 2
    save(strcat('simulation_model',num2str(model),'_expe',num2str(whichexpe),'_transfer'));
elseif phase_fit == 3
    save(strcat('simulation_model',num2str(model),'_expe',num2str(whichexpe),'_both'));
end


