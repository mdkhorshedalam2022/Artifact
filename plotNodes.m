function plotNodes(x,z)
%show distribution of nodes based on their types
global nodetype P
coordinates = [x(:) z(:)];
    for ii = 1:P
        plot(x(:),z(:),'o','MarkerSize',15,'Color','k');
        text(coordinates(ii,1),coordinates(ii,2),num2str(nodetype(ii)),'fontsize',5,'Color','r','fontweight','bold'); %node types
    end 
title('Distribution of Nodes','FontSize',20)
dim = [0.8 0.6 0.3 0.3];
str={'0: Core Nodes', '1: Bottom Left Node', '2: Bottom Right Node', '3: Top Left Node', '4: Top Right Node', '5: Bottom Nodes', '6: Top Nodes', '7: Left Side Nodes', '8: Right Side Nodes', '9: Inlet Nodes', '10: Oulet Nodes'};
annotation('textbox',dim,'String',str,'Color','r','FontSize',9,'FitBoxToText','on');
end