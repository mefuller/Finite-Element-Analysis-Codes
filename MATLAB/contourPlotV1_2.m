function contourPlotV1_2(comp, bottomlim, toplim, field, nodes, elements, snodes, enodes, Iorder, U, S, E)
% 20200629 - Added bottomlim and toplim inputs to supply to the colorbar
% maker function.  Now supplying the colorbar maker with mean element
% values instead of the entire array of nodal values.

    deformflag = 1;
    
    labelstr = {'x', 'y', 'xy'};
    
    if field == 'E'
        figure('Name', ['Strain ' labelstr{comp}]);
        title(sprintf(['Strain Contour Plot: \\epsilon_{' labelstr{comp} '}']));
    elseif field == 'S'
        figure('Name', ['Stress ' labelstr{comp}]);
        title(sprintf(['Stress Contour Plot: \\sigma_{' labelstr{comp} '}']));
    end

    xlabel('x', 'FontSize', 12); 
    ylabel('y', 'FontSize', 12);

    
    hold on
    for elem = 1:size(elements,2)
        n1=elements(1,elem);
        n2=elements(2,elem);
        n3=elements(3,elem);
        n4=elements(4,elem);
        x1 = nodes(1,n1) + deformflag*U(n1*2-1);
        x2 = nodes(1,n2) + deformflag*U(n2*2-1);
        x3 = nodes(1,n3) + deformflag*U(n3*2-1);
        x4 = nodes(1,n4) + deformflag*U(n4*2-1);
        y1 = nodes(2,n1) + deformflag*U(n1*2);
        y2 = nodes(2,n2) + deformflag*U(n2*2);
        y3 = nodes(2,n3) + deformflag*U(n3*2);
        y4 = nodes(2,n4) + deformflag*U(n4*2);
        if field == 'E'
            Se = E(:,elem*Iorder^2-(Iorder^2-1):elem*Iorder^2);
        elseif field == 'S'
            Se = S(:,elem*Iorder^2-(Iorder^2-1):elem*Iorder^2);
        end
        p = mean(Se(comp,:))*ones(2);
        pcolor([x4 x3;x1 x2], [y4 y3;y1 y2], p)
    end

    levels = 12;
    if field == 'E'
        cbar = abaquscolorbarV1_4(mean(enodes(comp,:,:),2),bottomlim,toplim,levels);
    elseif field == 'S'
        cbar = abaquscolorbarV1_4(mean(snodes(comp,:,:),2),bottomlim,toplim,levels);
    end
    axis equal
    hold off
end