clear

%% Load data

whichexpe = listdlg('SelectionMode','single','ListString',{'expe1','expe2','expe3'});

if whichexpe == 1
    load data_fig_expe1
elseif whichexpe == 2
    load data_fig_expe2
elseif whichexpe == 3
    load data_fig_expe3
end

colorbar = [125 0 250 80]/255;


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

%% Plot figure (opacity can vary with Matlab verion)

if whichexpe == 1

    order = 10:-1:1;

    figure;
    % Subjective values
    violinplotSB(explicit(order,:),Colors2(order,:),0,100);
    % Objective values
    rectangle('Position',[.5,13,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[1.5,49,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[2.5,13,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[3.5,31,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[4.5,49,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[5.5,13,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[6.5,85,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[7.5,13,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[8.5,49,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[9.5,85,1,2],'FaceColor',colorbar,'LineStyle','none');
    yline([100/3 50],'k:');
    set(gca,'Box','off') ;

    difference = {explicit(10,:)-14, explicit(8,:)-14, explicit(5,:)-14, explicit(3,:)-14,...
        explicit(7,:)-32,...
        explicit(9,:)-50, explicit(6,:)-50, explicit(2,:)-50,...
        explicit(4,:)-86, explicit(1,:)-86};

    figure;
    % Sub - Obj values
    violinplotSB(difference',Colors2([10 8 5 3 7 9 6 2 4 1],:),-80,60);
    yline(0,'k');
    set(gca,'Box','off');

elseif whichexpe == 2

    order = 10:-1:1;

    figure;
    % Subjective values
    violinplotSB(explicit(order,:),Colors2(order,:),0,100);
    % Objective values
    rectangle('Position',[.5,49,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[1.5,85,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[2.5,49,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[3.5,67,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[4.5,85,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[5.5,13,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[6.5,85,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[7.5,13,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[8.5,49,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[9.5,85,1,2],'FaceColor',colorbar,'LineStyle','none');
    yline([100/3 50],'k:');
    set(gca,'Box','off') ;

    difference_big = {explicit(5,:)-14, explicit(3,:)-14,...
        explicit(10,:)-50, explicit(8,:)-50, explicit(2,:)-50,...
        explicit(7,:)-68,...
        explicit(9,:)-86, explicit(6,:)-86, explicit(4,:)-86, explicit(1,:)-86};

    figure;
    % Sub - Obj values
    skylineplot_test(difference_big',Colors2 ([5 3 10 8 2 7 9 6 4 1],:),-80,60,12,'','','');
    set(gca,'Box','off');

else

    order = 12:-1:1;

    figure;
    skylineplot_test(explicit(order,:),Colors2(order,:),0,100,12,'','','','');
    rectangle('Position',[.5, 13,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[1.5,31,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[2.5,49,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[3.5,13,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[4.5,31,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[5.5,49,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[6.5,13,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[7.5,49,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[8.5,85,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[9.5,13,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[10.5,49,1,2],'FaceColor',colorbar,'LineStyle','none');
    rectangle('Position',[11.5,85,1,2],'FaceColor',colorbar,'LineStyle','none');
    yline([100/3 50],'k:');
    set(gca,'Box','off') ;

    difference = {explicit(12,:)-14, explicit(9,:)-14,explicit(6,:)-14, explicit(3,:)-14,...
        explicit(11,:)-32, explicit(8,:)-32,...
        explicit(10,:)-50, explicit(7,:)-50, explicit(5,:)-50, explicit(2,:)-50,...
        explicit(4,:)-86, explicit(1,:)-86};

    figure;
    skylineplot_test(difference',Colors2 ([12 9 6 3 11 8 10 7 5 2 4 1],:),-80,60,12,'','','');
    set(gca,'Box','off');

end


