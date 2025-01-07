clear


%%

load data_fig_expe3

subjects = 1:200;

rating_m=zeros(12,numel(subjects));


%% Simulate transfer Test

repet = 100;
POST = [unique(nchoosek(repmat(1:12,1,2),2),'rows');unique(nchoosek(repmat(1:12,1,2),2),'rows')]; POST(POST(:,1)==POST(:,2),:)=[];

for n = 1:repet
    
    clear Proba_post
    
    for sub = subjects
        
        ss = POST;
        Qf = explicit(:,sub);
        
        for j = 1:length(ss)
            
            % softmax
%             Proba_post{sub}(j)  = 1/(1+exp((Qf(ss(j,1)) - Qf(ss(j,2))) * .1));
            
            % argmax
            if Qf(ss(j,1)) ~= Qf(ss(j,2))
                [~, Proba_post{sub}(j)]  = max([Qf(ss(j,1)), Qf(ss(j,2))]);
            else
                Proba_post{sub}(j)  = .5;
            end
            
            % greedy
%             if rand < 0.95
%                 [~, Proba_post{sub}(j)]  = max([Qf(ss(j,1)), Qf(ss(j,2))]);
%             else
%                 Proba_post{sub}(j)  = .5;
%             end
            
        end
        
        Proba_post{sub} = Proba_post{sub}';
        
        % Transfer test
        
        Post_Test{1,sub} =[Proba_post{sub}(POST(:,2)==1,:)-1;(Proba_post{sub}(POST(:,1)==1,:)-2)*-1];
        Post_Test{2,sub} =[Proba_post{sub}(POST(:,2)==2,:)-1;(Proba_post{sub}(POST(:,1)==2,:)-2)*-1];
        Post_Test{3,sub} =[Proba_post{sub}(POST(:,2)==3,:)-1;(Proba_post{sub}(POST(:,1)==3,:)-2)*-1];
        Post_Test{4,sub} =[Proba_post{sub}(POST(:,2)==4,:)-1;(Proba_post{sub}(POST(:,1)==4,:)-2)*-1];
        Post_Test{5,sub} =[Proba_post{sub}(POST(:,2)==5,:)-1;(Proba_post{sub}(POST(:,1)==5,:)-2)*-1];
        Post_Test{6,sub} =[Proba_post{sub}(POST(:,2)==6,:)-1;(Proba_post{sub}(POST(:,1)==6,:)-2)*-1];
        Post_Test{7,sub} =[Proba_post{sub}(POST(:,2)==7,:)-1;(Proba_post{sub}(POST(:,1)==7,:)-2)*-1];
        Post_Test{8,sub} =[Proba_post{sub}(POST(:,2)==8,:)-1;(Proba_post{sub}(POST(:,1)==8,:)-2)*-1];
        Post_Test{9,sub} =[Proba_post{sub}(POST(:,2)==9,:)-1;(Proba_post{sub}(POST(:,1)==9,:)-2)*-1];
        Post_Test{10,sub}=[Proba_post{sub}(POST(:,2)==10,:)-1;(Proba_post{sub}(POST(:,1)==10,:)-2)*-1];
        Post_Test{11,sub}=[Proba_post{sub}(POST(:,2)==11,:)-1;(Proba_post{sub}(POST(:,1)==11,:)-2)*-1];
        Post_Test{12,sub}=[Proba_post{sub}(POST(:,2)==12,:)-1;(Proba_post{sub}(POST(:,1)==12,:)-2)*-1];
        
        for m=1:12
            rating_m(m,sub)=rating_m(m,sub)+nanmean(Post_Test{m,sub})/repet;
        end
    end
end


%%


order = [12 11 10 9 8 7 6 5 4 3 2 1];


%%


figure;
bars_datamodel(rating_m(order,subjects),rating(order,subjects),Colors2(order,:),0,1,8,'','','');
ylim([0 1]);
yline(0.5,'k:');
yline(1/3,'k:');
set(gca,'Box','off');



%%

figure;
scatterCorrSpear(mean(rating(:,subjects),2),mean(rating_m(:,subjects),2),1);


%%

figure;
for i = order
    errorbar(mean(rating(i,subjects)),mean(rating_m(i,subjects)),sem(rating(i,subjects)),sem(rating_m(i,subjects)),"both","o",...
        'MarkerFaceColor',Colors2(i,:),'MarkerEdgeColor',Colors2(i,:),'Color',Colors2(i,:));
    hold on
end
xlim([0 1]);
ylim([0 1]);
line([0 1],[0 1],'Color','k');
pbaspect([1 1 1]);
