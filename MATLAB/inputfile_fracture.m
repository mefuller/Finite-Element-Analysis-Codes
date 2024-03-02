jobfolder = pwd;

%% ########## Material constants and orientation ##########
E_m = 1e5;
nu_m = 0.1;
Ax = 1e5;
Bx = -1e5;
theta0 = 2/3*pi;
E_f = 1e5;
nu_f = 0.1;
v_f = 0.1;
mat_ndir = 2;

%% ########## Load setup ###############
trac = 5e3; trac_dir = 2;

%% ########## Analysis controls ###############
maxit =  10; tol = 1e-3; nsteps = (50)*trac; Iorder = 2;

%% ########## Output Setup ##################
intervals = 60;
movieframerate = 12;
histfreq = floor(nsteps/intervals); %history request intervals
frames = 1:histfreq:nsteps;
if ~ismember(nsteps,frames)
    frames = [frames nsteps];
end

domflag = 0; %0 or 1; Makes a plot of the domain with node and gauss point numbers.

%% ############ Mesh Setup ################
ndivs = 6;
numrings = 40;
numrect = 10;
origin = [1 0];
xb = [0 2];
yb = [0 1 5];
[nodes, elements, adj, tracsetup, rhos] = makeMesh(xb, yb, numrings, ndivs, numrect, origin);

%% ############ Boundary Conditions Setup ################
%Boundary Conditions.  DOFs is a list of DOFs that are to be constrained.
% A is the prescribed values at those DOFs:
DOFs = [find(nodes(1,:) >= 1 & nodes(2,:) == 0)*2, ...
    find(nodes(1,:) == 1 & nodes(2,:) == 0)*1];   %DOFs involved in BC
A =    (zeros(1,length(DOFs)));   %Prescribed values at DOFs

% nodes(1,:) = nodes(1,:) - 1;
% origin = [0 0];

% ndivs = [16 16];
% xb = [0 2];
% yb = [0 2];
% [nodes, elements, adj, tracsetup] = makeRectMesh(xb, yb, ndivs);




% Test mesh
% DOFs = [1 2 4 7];
% A = [0 0 0 0];
% nodes = [0 0; 1 0; 1 1; 0 1;]';
% elements = [1 2 3 4]';
% tracsetup = [3 4];