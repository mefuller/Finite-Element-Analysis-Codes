function [ snodes ] = extrapnodesV1_4(Iorder, numnodes, Smatrix, elements )
%Extrapolates stresses to nodes.
%   Detailed explanation goes here
    %Extrapolation matrix for gauss points to nodes
    elements = elements;
    S = Smatrix;
    
    if Iorder == 1;
        Se = [];
        snodes = [];
        for elem = 1:size(elements,2);  %Perform extrapolation for all elements
            Se(1:4,elem) = S(1,elem)*ones(4,1); %s11
            Se(5:8,elem) = S(2,elem)*ones(4,1); %s22
            Se(9:12,elem) = S(3,elem)*ones(4,1); %s12
            snodes(1,:,elem) = Se(1:4,elem)';
            snodes(2,:,elem) = Se(5:8,elem)';
            snodes(3,:,elem) = Se(9:12,elem)';
        end
    end
    if Iorder == 2

        extrap = [1+sqrt(3)/2 -1/2 1-sqrt(3)/2 -1/2;...
                  -1/2 1+sqrt(3)/2 -1/2 1-sqrt(3)/2;...
                  1-sqrt(3)/2 -1/2 1+sqrt(3)/2 -1/2;...
                  -1/2 1-sqrt(3)/2 -1/2 1+sqrt(3)/2];

        %perform of stresses extrapolation to nodes      
        Se = [];
        snodes = [];
        for elem = 1:size(elements,2);  %Perform extrapolation for all elements
            Se(1:4,elem) = (extrap*S(1,elem*4-3:elem*4)')'; %s11
            Se(5:8,elem) = (extrap*S(2,elem*4-3:elem*4)')'; %s22
            Se(9:12,elem) = (extrap*S(3,elem*4-3:elem*4)')'; %s12
            snodes(1,:,elem) = Se(1:4,elem)';
            snodes(2,:,elem) = Se(5:8,elem)';
            snodes(3,:,elem) = Se(9:12,elem)';
        end
    end
    
end

