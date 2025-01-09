clear

%% SET UP

% COLUMN 1: participant number
% COLUMN 2: phase number (0 training, 1 learning, 2 transfer)
% COLUMN 3: trial number
% COLUMN 4: context number
%     expe 1
%       1 = WT (86/50/14)
%       2 = WB (86/-/14)
%       3 = NT (50/32/14)
%       4 = NB (50/-/14)
%     expe 2
%       1 = WT (86/50/14)
%       2 = WB (86/-/14)
%       3 = NT (86/68/50)
%       4 = NB (86/-/50)
%     expe 3
%       1 = WT (86/50/14)
%       2 = WT (86/50/14)
%       3 = NT (50/32/14)
%       4 = NT (50/32/14)
% COLUMN 5: left option
% COLUMN 6: middle option
% COLUMN 7: right option
% COLUMN 8: choice -1 left, 0 middle, 1right
% COLUMN 9: choice accuracy
% COLUMN 10: outcome chosen option
% COLUMN 11: outcome unchosen option 1
% COLUMN 12: outcome unchosen option 2
% COLUMN 13: trial reaction time (ms)
% COLUMN 14: participant ID

whichexpe = listdlg('SelectionMode','single','ListString',{'expe1','expe2','expe3'});

if whichexpe == 1
    load data_expe1
    load data_expe1_explicit
elseif whichexpe == 2
    load data_expe2
    load data_expe2_explicit
elseif whichexpe == 3
    load data_expe3
    load data_expe3_explicit
end

subjects = 1:max(data(:,1));

sub=0;
for nsub = subjects
    
    sub=sub+1;
    
    % Accuracy per condition
    if whichexpe < 3
        cond1(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)==1,9);
        cond2(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)==2,9);
        cond3(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)==3,9);
        cond4(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)==4,9);
    else
        cond1(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)==1,9);
        cond2(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & abs(data(:,4))==2,9);
        cond3(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)==3,9);
        cond4(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & abs(data(:,4))==4,9);
    end
    % 
    % % Overall accuracy learning / transfer
    % perfLT(:,sub) = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)>0,9);
    % perfTT(:,sub) = data(data(:,1)==nsub & data(:,2)==2,9);
    % 
    % % Learning options and choices
    % symLT{sub}  = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)>0,5:7);
    % choLT{sub}  = data(data(:,1)==nsub & data(:,2)==1 & data(:,4)>0,8)+2;

    % Trial options
    sym{sub}  = data(data(:,1)==nsub & data(:,2)==2,5:7);
    % Trial choice
    cho{sub}  = data(data(:,1)==nsub & data(:,2)==2,8);

    % For each option, choice rate and explicit rating
    for i=1:(10+2*(whichexpe==3))
        rating(i,sub) = (5+1*(whichexpe==3)) * nanmean((sym{sub}(:,1)==i& cho{sub}==-1) + (sym{sub}(:,2)==i& cho{sub}==0) + (sym{sub}(:,3)==i& cho{sub}==1));
        explicit(i,sub) = nanmean(dataE(dataE(:,1)==sub & dataE(:,4)==i-1,5));
    end
    
    % Check for incomplete datasets (different trial number for expe3)
    if length(data(data(:,1)==nsub & data(:,2)==2,:)) ~= (180 - 48*(whichexpe==3))
        disp(sub);
    end

    % Check for incomplete datasets (different trial number for expe3)
    if length(dataE(dataE(:,1)==sub,:)) ~= (40 - 16*(whichexpe==3))
        disp('missing participants in explicit phase'); disp(sub);
    end
    
end

%% Save data for figures

if whichexpe == 1
    save data_fig_expe1
elseif whichexpe == 2
    save data_fig_expe2
elseif whichexpe == 3
    save data_fig_expe3
end
