%% Nonlinear FEM code v1_8
% Plane Strain Formulation
% by N. Payne 2016 - 2020

%This function runs the inputfile script designated by the string 'inp' to
%initialize the analysis inputs and then returns U E S and S33
%for the FEA problem.

%% Notes
% 1.) S33 introduced to collect stress component in the z direction.
% Currently not used.
% 2.) Introduced Runge Kutta 4th Order used to compute stress increment.

function [U,E,S,S33,history] = sol_nlFEM_18(inp)
%% Input of Parameters

rev_date = '2-Aug-2020';
disp(['timedate:' string(datetime('now'))])
disp(['**Analysis information**' newline...
    'Nonlinear FEM, Plane strain formulation.' newline...
    'Solver type:  Quasi-Newton algorithm' newline...
    'sol_nlFEM_18() solver version 1.8; rev. date: ' rev_date])



eval(inp)
logfile = fopen('logfile.txt','w');
warning('off','MATLAB:nearlySingularMatrix')

history = cell(3,1);


%% Gauss Point Generation
GP{1} = [0 0 2*2];

GP{2} = [-1/sqrt(3) -1/sqrt(3) 1;
      1/sqrt(3) -1/sqrt(3) 1;
      1/sqrt(3) 1/sqrt(3) 1;
      -1/sqrt(3) 1/sqrt(3) 1];
  
GP{3} = [-sqrt(3/5) -sqrt(3/5) (5/9)*(5/9);...
         0 -sqrt(3/5) (8/9)*(5/9);...
         sqrt(3/5) -sqrt(3/5) (5/9)*(5/9);...
         -sqrt(3/5) 0 (5/9)*(8/9);...
         0 0 (8/9)*(8/9);...
         sqrt(3/5) 0 (5/9)*(8/9);...
         -sqrt(3/5) sqrt(3/5) (5/9)*(5/9);...
         0 sqrt(3/5) (8/9)*(5/9);...
         sqrt(3/5) sqrt(3/5) (5/9)*(5/9);];

GPweight = [2 1 8/9 5/9];

for i = 1:Iorder^2
    xi = GP{Iorder}(i,1); eta = GP{Iorder}(i,2);
    N(i*2-1:i*2,:) = [(1/4)*(1-xi).*(1-eta) 0  (1/4)*(1+xi).*(1-eta) 0  (1/4)*(1+xi).*(1+eta) 0  (1/4)*(1-xi).*(1+eta) 0;...
         0  (1/4)*(1-xi).*(1-eta) 0  (1/4)*(1+xi).*(1-eta) 0  (1/4)*(1+xi).*(1+eta) 0  (1/4)*(1-xi).*(1+eta)];
end
 
nume = size(elements,2);
Q = nodes(:,elements(:)); Q = reshape(Q,[8,nume]);
gps = N*Q; gps = reshape(gps, [2,nume*Iorder^2]);
gps(3,:) = repmat(GP{Iorder}(:,3)',[1 size(elements,2)]);

%%
%Calculate B and detJ for each gauss point.
for elem = 1:size(elements,2)
    for i = 1:Iorder^2
        [B(:,:,elem*Iorder^2-(Iorder^2-i)), detJ(elem*Iorder^2-(Iorder^2-i))]= BmatdetJ(reshape(nodes(:,elements(:,elem)),[8 1]),GP{Iorder}(i,:));
    end
end
% reshape(nodes(:,elements(:,2)),[8 1])
%% Initialize Matrices

%Displacement and displacement increment:
U = zeros(length(nodes)*2,1);
delA = zeros(length(nodes)*2,1);

%Stress, stress increment, strain, strain increment:
S = zeros(3,size(gps,2)); 
E = zeros(3,size(gps,2));
nS = zeros(3,size(gps,2));
S33 = zeros(1,size(gps,2));
nE = zeros(3,size(gps,2)); 

%Internal Force and Traction Vector and Global Stiffness Matrix
T = sparse(length(nodes)*2,1);
F = sparse(length(nodes)*2,1);
K = sparse(length(nodes)*2,length(nodes)*2);

%Point force boundary conditions
%Set traction boundary condition on upper face of domain.


for k = 1:size(tracsetup,1)
    n1 = tracsetup(k,1);
    x1 = nodes(1,n1);
    y1 = nodes(2,n1);
    n2 = tracsetup(k,2);
    x2 = nodes(1,n2);
    y2 = nodes(2,n2);
    len23 = sqrt((x2-x1)^2 + (y2-y1)^2);
    ty = 0;
    tx = 0;
    if trac_dir == 1
        tx = trac;
    elseif trac_dir == 2
        ty = trac;
    end
    F(n1*2-1) = F(n1*2-1) + tx*len23/2;
    F(n1*2) = F(n1*2) + ty*len23/2;
    F(n2*2-1) = F(n2*2-1) + tx*len23/2;
    F(n2*2) = F(n2*2) + ty*len23/2;
end


%% Plot of Meshed Domain

if domflag == 1
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

    %Showing gauss points for just first element:
    figure('Name', 'FE Domain')
    gplot(adj, nodes')
    hold on
    
    %scatter(gps(1,:),gps(2,:))
    title(sprintf('FE Domain and Node Numbers'));
    xlabel('x', 'FontSize', 12); 
    ylabel('y', 'FontSize', 12);

    for i = 1:size(gps,2)
        text(gps(1,i), gps(2,i), num2str(i),...
            'Color', 'red')
    end

    for i = 1:length(nodes)
        text(nodes(1,i), nodes(2,i), num2str(i),...
            'Color', 'blue')
    end
    hold off
end
%% Setup of Local to Global Matrix Mappings 
Ig = zeros(64*size(elements,2),1);
Jg = zeros(64*size(elements,2),1);
IgT = zeros(8*size(elements,2),1);
JgT = zeros(8*size(elements,2),1);
cc = 1;

for elem = 1:size(elements,2)
    mapping = [1 2 3 4 5 6 7 8;...
    elements(1,elem)*2-1 elements(1,elem)*2 ...
    elements(2,elem)*2-1 elements(2,elem)*2 ...
    elements(3,elem)*2-1 elements(3,elem)*2 ...
    elements(4,elem)*2-1 elements(4,elem)*2];
    [II, JJ] = meshgrid(mapping(2,:),mapping(2,:));
    Ig(cc*64-63:cc*64) = II(:);
    Jg(cc*64-63:cc*64) = JJ(:);
    IgT(cc*8-7:cc*8) = mapping(2,:);
    JgT(cc*8-7:cc*8) = ones(8,1);
    cc = cc + 1;
end
%% Solve Block: Newton-Raphson Method
tic; fprintf(logfile,'Calculations started.');
for step = 1:nsteps
    fprintf(logfile,'Step:  ',num2str(step));
    error = 1.0; nit = 0; factor = step/nsteps; errflag = 0;
    
    delA = zeros(length(nodes)*2,1);
    delE = zeros(3,size(gps,2));
    delS = zeros(3,size(gps,2));
       
    while nit < maxit
            nit = nit + 1; fprintf(logfile,['\n' 'step:  ' num2str(step) '; iteration:  ' num2str(nit)]);
            if nit == 1
                K = sparse(length(nodes)*2,length(nodes)*2);
                Kg = zeros(64*size(elements,2),1);
            end
            cc = 1;
            T = sparse(length(nodes)*2,1);
            Tg = zeros(8*size(elements,2),1);
            for elem = 1:size(elements,2)
                k_e = zeros(8);
                T_e = zeros(8,1);
                n1=elements(1,elem);
                n2=elements(2,elem);
                n3=elements(3,elem);
                n4=elements(4,elem);
                delU_e = [delA(n1*2-1:n1*2); delA(n2*2-1:n2*2);...
                    delA(n3*2-1:n3*2); delA(n4*2-1:n4*2)];
                for gp = elem*Iorder^2-(Iorder^2-1):elem*Iorder^2
                    Cm(:,:,gp) = Cmat_transvIsoV3b(mat_ndir, E(:,gp), E_m, nu_m, Ax, Bx, theta0, E_f, nu_f, v_f);
                    weight = gps(3,gp);
                    Bm = B(:,:,gp);
                    dJ = detJ(gp);
                    
                    delE(:,gp) = Bm*delU_e;
                    nE(:,gp) = E(:,gp) + delE(:,gp);
                    
                    nS(:,gp) = S(:,gp) + RK4V2_delS(delE(:,gp), E(:,gp), mat_ndir, E_m, nu_m, Ax, Bx, theta0, E_f, nu_f, v_f);
                    
                    T_e = T_e + transpose(Bm)*(nS(:,gp))*weight*dJ;
                    if nit == 1
                        k_e = k_e + transpose(Bm)*Cm(:,:,gp)*Bm*weight*dJ;
                    end
                end
                %Scattering Ke to Kglobal and Te to Tglobal:
                mapping = [1 2 3 4 5 6 7 8;...
                    elements(1,elem)*2-1 elements(1,elem)*2 ...
                    elements(2,elem)*2-1 elements(2,elem)*2 ...
                    elements(3,elem)*2-1 elements(3,elem)*2 ...
                    elements(4,elem)*2-1 elements(4,elem)*2];
                if nit == 1
                    Kg(cc*64-63:cc*64) = k_e(:);                         
                end
                Tg(cc*8-7:cc*8) = T_e(:);
                cc = cc + 1;
            end
            if nit == 1
                K = sparse(Ig,Jg,Kg);
            end
            T = sparse(IgT,JgT,Tg);


        %Compute Residual applied nodal force BC - internal nodal forces.
%         factor*F
%         T
        R = factor*F - T;
%% Boundary Conditions 
        
%         DOFs = [find(nodes(1,:) == 0)*2-1, ...
%             find(nodes(2,:) == 0)*2];   %DOFs involved in BC
%         A =    (zeros(1,length(DOFs)));   %Prescribed values at DOFs
        
       %% Elimination Approach 
       Krows = [];  %Matrix to hold K terms of BC DOFs
       activeDOFs = 1:length(K);
       activeDOFs(DOFs)=[];   %Active DOFs in model after BCs applied.
              
       for i = 1:length(DOFs)
          p = DOFs(i);
          Krows(i,:) = K(p,:);  %Set K rows aside for later
          
          K(p,:) = 0; %delete p row
          K(:,p) = 0; %delete p column
          K(p,p) = 1;
          R(p) = A(i);
       end
       for k = 1:length(activeDOFs)
          i = activeDOFs(k);
          modF = 0;
          for r = 1:length(DOFs)
              modF = modF + Krows(r,i)*A(r);
          end
          R(i) = R(i) - modF;   %Modify Residual vector to enforce BC
       end
        
        %% Lagrange Multiplier Approach
%         C = zeros(length(DOFs),length(nodes)*2);
%         Q = [];
%         for i = 1:length(DOFs)
%             C(i,DOFs(i)) = 1;
%             Q(i,1) = A(i);
%         end       
       %%
        %Check for convergence: err < tol;
%         error = sqrt(R'*R)/sqrt(factor*F'*factor*F);
%         error = sqrt(R'*R)/sqrt(F'*F);
        error = max(R)/mean(F);
        fprintf(logfile,['Rerr: ' num2str(error) ' ']);
        if error < tol
            errflag = 1;
        else
            errflag = 0;
        end
        
        if errflag == 1
            break
        end
        
        da = K\R;
       
%         dadlamb = [K C';C zeros(size(C,1))]\[R;Q];
%         da = dadlamb(1:length(nodes)*2);
%         lamb = dadlamb(length(nodes)*2+1:length(dadlamb));
        
        delA = delA + da;
        
        %Displacment norm err indicator.
%         err = sqrt((transpose(da)*da)/(transpose(U+delA)*(U+delA)))
%             error = max(da)/max(U);
%             fprintf(logfile,['Cerr: ' num2str(error) ' ']);
% 
%             if error < tol
%                 errflag = 1;
%             else
%             errflag = 0;
%             end
            
    end
    
	if nit == maxit
        error('Max. iterations reached')
    end
    
 	if imag(error) ~= 0
        error('Imaginary number in solution.')
    end
    
    %Update strains, stresses, and displacements after convergence.
    E = nE;
    S = nS;
    U = U + delA;
    
    %Write history outputs:
    if ismember(step,frames)
        history{1,step} = E;
        history{2,step} = S;
        history{3,step} = U;
    end
    
    %Display progress percentage:
    if ismember(step,0:nsteps/20:nsteps)
        disp([num2str(step/nsteps*100) '%'])
    end
end
toc; fprintf(logfile,'Calculations completed.');
fclose(logfile);
warning('on','MATLAB:nearlySingularMatrix');
end 