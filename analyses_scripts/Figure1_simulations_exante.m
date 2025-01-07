clear


%% SET UP

model1 = 1; % Q-LEARNING
model2 = 2; % DIVISIVE NORMALIZATION
model3 = 3; % RANGE NORMALIZATION

whichexpe = listdlg('SelectionMode','single','ListString',{'expe1','expe2'});

if whichexpe == 1
    sim1 = load(strcat('exante_model',num2str(model1),'_expe1'));
    sim2 = load(strcat('exante_model',num2str(model2),'_expe1'));
    sim3 = load(strcat('exante_model',num2str(model3),'_expe1'));
else
    sim1 = load(strcat('exante_model',num2str(model1),'_expe2'));
    sim2 = load(strcat('exante_model',num2str(model2),'_expe2'));
    sim3 = load(strcat('exante_model',num2str(model3),'_expe2'));
end

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


%%

order = [10 9 8 7 6 5 4 3 2 1];
Colors2 = Colors2(order,:);

simu = [mean(sim1.rating_arg(order,:),2)'; mean(sim2.rating_arg(order,:),2)'; mean(sim3.rating_arg(order,:),2)'];


%%

figure;

b = bar(simu, 'LineWidth',1.5,'FaceAlpha',0.7);

hold on

for i=1:numel(b)

    xData = b(i).XData+b(i).XOffset;

    % bar color
    b(i).FaceColor = Colors2(i,:);
    b(i).EdgeColor = Colors2(i,:);
    hold on

end

yline([0.5 1/3],'k:')
ylim([0 1])

pbaspect([3 1 1]);
set(gcf,'Position',[100 200 1300 400])
set(gca, 'Box', 'off') ;











