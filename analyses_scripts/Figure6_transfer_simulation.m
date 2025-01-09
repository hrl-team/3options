clear


%%

whichexpe = listdlg('SelectionMode','single','ListString',{'expe1','expe2','expe3'});

if whichexpe == 1
    load data_fig_expe1
    rating_m=zeros(10,numel(subjects));
elseif whichexpe == 2
    load data_fig_expe2
    rating_m=zeros(10,numel(subjects));
elseif whichexpe == 3
    load data_fig_expe3
    rating_m=zeros(12,numel(subjects));
end


%% Choose colors for figures

if whichexpe < 3

    Colors(1,:) = [135 62 35] /255;
    Colors(2,:) = [21 76 121] /255;
    Colors(3,:) = [226 135 67]/255;
    Colors(4,:) = [30 129 176] /255;

    Colors2(1,:) = Colors(1,:);
    Colors2(2,:) = Colors(1,:);
    Colors2(3,:) = Colors(1,:);
    Colors2(4,:) = Colors(2,:);
    Colors2(5,:) = Colors(2,:);
    Colors2(6,:) = Colors(3,:);
    Colors2(7,:) = Colors(3,:);
    Colors2(8,:) = Colors(3,:);
    Colors2(9,:) = Colors(4,:);
    Colors2(10,:) = Colors(4,:);

else

    Colors(1,:) = [29 81 90]  /255;
    Colors(2,:) = [97 179 159] /255;
    Colors(3,:) = [159 198 148] /255;
    Colors(4,:) = [224 133 89] /255;

    Colors2(1,:) = Colors(1,:);
    Colors2(2,:) = Colors(1,:);
    Colors2(3,:) = Colors(1,:);
    Colors2(4,:) = Colors(2,:);
    Colors2(5,:) = Colors(2,:);
    Colors2(6,:) = Colors(2,:);
    Colors2(7,:) = Colors(3,:);
    Colors2(8,:) = Colors(3,:);
    Colors2(9,:) = Colors(3,:);
    Colors2(10,:) = Colors(4,:);
    Colors2(11,:) = Colors(4,:);
    Colors2(12,:) = Colors(4,:);

end

%% Simulate transfer Test

% Number of iterations
repet=100;

% Create the transfer 180x2 or 132x2 matrix (all combinations)
if whichexpe == 3
    POST = [unique(nchoosek(repmat(1:10,1,2),2),'rows');unique(nchoosek(repmat(1:10,1,2),2),'rows')]; POST(POST(:,1)==POST(:,2),:)=[];
else
    POST = [unique(nchoosek(repmat(1:12,1,1),2),'rows');unique(nchoosek(repmat(1:12,1,1),2),'rows')]; POST(POST(:,1)==POST(:,2),:)=[];
end

for n = 1:repet

    clear simchoTT

    for sub = subjects

        post_sym = POST;
        Qfinal   = explicit(:,sub);

        for j = 1:length(post_sym)

            % softmax
            %             simchoTT{sub}(j)  = 1/(1+exp((Qf(ss(j,1)) - Qf(ss(j,2))) * .1));

            % argmax
            if Qfinal(post_sym(j,1)) ~= Qfinal(post_sym(j,2))
                [~, simchoTT{sub}(j)]  = max([Qfinal(post_sym(j,1)), Qfinal(post_sym(j,2))]);
            else
                simchoTT{sub}(j)  = .5;
            end

            % greedy
            %             if rand < 0.95
            %                 [~, simchoTT{sub}(j)]  = max([Qf(ss(j,1)), Qf(ss(j,2))]);
            %             else
            %                 simchoTT{sub}(j)  = .5;
            %             end

        end

        simchoTT{sub} = simchoTT{sub}';

        % Transfer test

        choice_TT{1,sub} =[simchoTT{sub}(POST(:,2)==1,:)-1;(simchoTT{sub}(POST(:,1)==1,:)-2)*-1];
        choice_TT{2,sub} =[simchoTT{sub}(POST(:,2)==2,:)-1;(simchoTT{sub}(POST(:,1)==2,:)-2)*-1];
        choice_TT{3,sub} =[simchoTT{sub}(POST(:,2)==3,:)-1;(simchoTT{sub}(POST(:,1)==3,:)-2)*-1];
        choice_TT{4,sub} =[simchoTT{sub}(POST(:,2)==4,:)-1;(simchoTT{sub}(POST(:,1)==4,:)-2)*-1];
        choice_TT{5,sub} =[simchoTT{sub}(POST(:,2)==5,:)-1;(simchoTT{sub}(POST(:,1)==5,:)-2)*-1];
        choice_TT{6,sub} =[simchoTT{sub}(POST(:,2)==6,:)-1;(simchoTT{sub}(POST(:,1)==6,:)-2)*-1];
        choice_TT{7,sub} =[simchoTT{sub}(POST(:,2)==7,:)-1;(simchoTT{sub}(POST(:,1)==7,:)-2)*-1];
        choice_TT{8,sub} =[simchoTT{sub}(POST(:,2)==8,:)-1;(simchoTT{sub}(POST(:,1)==8,:)-2)*-1];
        choice_TT{9,sub} =[simchoTT{sub}(POST(:,2)==9,:)-1;(simchoTT{sub}(POST(:,1)==9,:)-2)*-1];
        choice_TT{10,sub}=[simchoTT{sub}(POST(:,2)==10,:)-1;(simchoTT{sub}(POST(:,1)==10,:)-2)*-1];
        if whichexpe == 3
            choice_TT{11,sub}=[simchoTT{sub}(POST(:,2)==11,:)-1;(simchoTT{sub}(POST(:,1)==11,:)-2)*-1];
            choice_TT{12,sub}=[simchoTT{sub}(POST(:,2)==12,:)-1;(simchoTT{sub}(POST(:,1)==12,:)-2)*-1];
        end

        for m=1:(10+2*(whichexpe==3))
            rating_m(m,sub)=rating_m(m,sub)+nanmean(choice_TT{m,sub})/repet;
        end

    end
end


%% Bar order

if whichexpe < 3
    order = 10:-1:1;
else
    order = 12:-1:1;
end


%% Plot figure

figure;
bars_datamodel(rating_m(order,:),rating(order,:),Colors2(order,:),0,1,8,'','','');
ylim([0 1]);
yline(0.5,'k:');
yline(1/3,'k:');
set(gca,'Box','off');


%% Supp. figure

figure;
for i = order
    errorbar(mean(rating(i,:)),mean(rating_m(i,:)),sem(rating(i,:)),sem(rating_m(i,:)),"both","o",...
        'MarkerFaceColor',Colors2(i,:),'MarkerEdgeColor',Colors2(i,:),'Color',Colors2(i,:));
    hold on
end
xlim([0 1]);
ylim([0 1]);
line([0 1],[0 1],'Color','k');
pbaspect([1 1 1]);
