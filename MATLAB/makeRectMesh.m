function [nodes, elements, adj, tracsetup] = makeRectMesh(xb, yb, divs);


    divs = divs; 
    origin = origin;
    xb = xb;
    yb = yb;  

    nodes = [];
    elements = [];


       %% Create rectangular portion
    %First, make transition layer
    for j = (yb(1):(yb(2)-yb(1))/divs(2):yb(2))
        nodes = [nodes; [(xb(1):(xb(2)-xb(1))/divs(1):xb(2))' repmat(j,[divs(1)+1,1])]];
    end
    for j = 0:divs(2)-1    
        for i = 1:divs(1)
            elements = [elements; i+j*(divs(1)+1) i+j*(divs(1)+1)+1 i+(j+1)*(divs(1)+1)+1  i+(j+1)*(divs(1)+1)];
        end
    end
    %% Plot of domain
    %Make adjacency matrix for visualization
    adj = sparse(length(nodes),length(nodes));
    for row = 1:size(elements,1);
        for col = 1:size(elements,2);
            i = elements(row,col);
            if col == size(elements,2);
                j = elements(row,1);
            else
                j = elements(row,col+1);
            end
            adj(i,j) = adj(i,j) + 1;
            adj(j,i) = adj(j,i) + 1;
        end
    end
    adj(find(adj > 1)) = 1;
    
    nodes = nodes';
    elements = elements';
    
    %% Selection surface for analysis loads
    edgenodes = find(nodes(2,:) == yb(2));
    edgeelements = [];
    for xe = edgenodes
        [row, col] = find(elements == xe);
        edgeelements = [edgeelements col'];
    end
    edgeelements = unique(edgeelements);
%     [row, col] = find(elements(1,:) == edgenodes' | elements(2,:) == edgenodes' | elements(3,:) == edgenodes' | elements(4,:) == edgenodes');
    tracsetup = [];
    for i = 1:length(edgeelements)
        for en1 = 1:length(edgenodes)
            for en2 = 1:length(edgenodes)
                if en1 == en2
                    break
                end
                if ismember(edgenodes(en1),elements(:,edgeelements(i))) == 1 & ismember(edgenodes(en2),elements(:,edgeelements(i))) == 1
                    tracsetup = [tracsetup; [edgenodes(en1) edgenodes(en2)]];
                end
               
            end
        end
    end

   