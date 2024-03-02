function [nodes, elements, adj, tracsetup, rhos] = makeMesh(xb, yb, numrings, divs, numrect, origin);


    %Sets to track
    rhos = [];
    thetas = [];

    numrings = numrings; %in crack section
    divs = divs; %per 45deg
    numrect = numrect;   %num elements in rectangular section
    origin = origin;
    xb = xb;
    yb = yb; %y0 crack, y1 crack, y2 rectangle   

    nodes = [origin];
    elements = [];

    Ri = 9E-9;

    Ntotal = numrings;
    i = 1:Ntotal-1;

    %Set rho log spacing
    rho_spacing = Ri*(10.^((8/(Ntotal-2))*(i-1)));

    %set angular spacing
    theta_inc = pi/4/divs;
    ringarray = [cos(0:theta_inc:pi)' sin(0:theta_inc:pi)'];
    ringarray(divs*2+1,1) = 0; ringarray(divs*4+1,2) = 0;


    %% Build first ring
    nodes = [nodes; (rho_spacing(1)*ringarray) + origin];
    
    
    rhos = 2:size(nodes,1); %add first ring
    
    for i = 1:divs*4
        elements = [elements; 1 i+1 i+2 1];
    end
    %% Build all other rings
    for rnum = 2:length(rho_spacing)
        l1 = length(nodes);
        nodes = [nodes; (rho_spacing(rnum)*ringarray) + origin];
        
        %add to rhos
        newnodes = l1 + 1:l1 + length(ringarray);
        rhos = [rhos; newnodes];
        
        nring = (divs*4+1);
        k = 1 + nring*(rnum-2);
        for i = 1:divs*4
            elements = [elements; [k+i k+i+nring k+i+nring+1 k+i+1]];
        end
    end
    
    %% Build outer square ring
    for rnum = length(rho_spacing)+1
        nodes = [nodes; (ringarray) + origin];
        nring = (divs*4+1);
        k = 1 + nring*(rnum-2);
        for i = 1:divs*4
            elements = [elements; [k+i k+i+nring k+i+nring+1 k+i+1]];
        end
    end
    %Move outer shell nodes
    edgenodes = (length(nodes)-(divs*4)):length(nodes);
    for i = 1:nring
        if i <= divs
            nodes(edgenodes(i),1) = xb(2);
            nodes(edgenodes(i),2) = (xb(2)-origin(1))*tan(theta_inc*(i-1));
        elseif divs < i && i < divs*3+1
            nodes(edgenodes(i),1) = (yb(2)-origin(2))/tan(theta_inc*(i-1))+origin(1);
            nodes(edgenodes(i),2) = yb(2);
        else
            nodes(edgenodes(i),1) = xb(1);
            nodes(edgenodes(i),2) = -(xb(2)-origin(1))*tan(theta_inc*(i-1));
        end
    end

    %% Create rectangular portion
    %First, make transition layer
    elemspace = (yb(3)-yb(2))/numrect;
    k = length(ringarray)-divs;

    i = length(ringarray)*numrings+1-(divs*3):length(ringarray)*numrings+1-divs;
    nodes = [nodes; [nodes(i,1) nodes(i,2) + elemspace]];

    i = length(ringarray)*numrings+1-(divs*3):length(ringarray)*numrings-divs;
    elements = [elements; [(i)'+1 i' i'+k i'+k+1]];

    %Now, make the rest of the rectangular section
    for p = 2:numrect
        i = length(ringarray)*numrings+1-(divs*3):length(ringarray)*numrings+1-divs;
        nodes = [nodes; [nodes(i,1) nodes(i,2) + elemspace*p]];
    end

    k = divs*2+1;
    for p = 0:numrect-2
        i = length(ringarray)*numrings+2+p*(divs*2+1):length(ringarray)*numrings+2*divs+1+p*(divs*2+1);
        elements = [elements; [(i)'+1 i' i'+k i'+k+1]];
    end

    nodes(find(nodes(:)<1e-12))=0; %clean small numbers up

    %% Plot of domain
    %Make adjacency matrix for visualization

    nodes = nodes';
    elements = elements';   
    
    adj = zeros(length(nodes));
    for elem = 1:size(elements,2);
        for node = 1:size(elements,1);
            i = elements(node,elem);
            if node == size(elements,1);
                j = elements(1,elem);
            else
                j = elements(node+1,elem);
            end
            adj(i,j) = adj(i,j) + 1;
            adj(j,i) = adj(j,i) + 1;
        end
    end
    adj(find(adj > 1)) = 1;

    
    %% Selection surface for analysis loads
    
    edgenodes = find(nodes(2,:) == yb(3));
    tracsetup = [];
    for i = 1:length(elements)
        for en1 = 1:length(edgenodes)
            for en2 = 1:length(edgenodes)
                if en1 == en2
                    break
                end
                if ismember(edgenodes(en1),elements(:,i)) == 1 && ismember(edgenodes(en2),elements(:,i)) == 1
                    tracsetup = [tracsetup; [edgenodes(en1) edgenodes(en2)]];
                end
            end
        end
    end
    
%     rhos = flip(rhos,2); %flip order so that theta = 0 is the column 1.