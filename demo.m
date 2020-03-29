% 
% Very short description of the files in this folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% mixturefit.m 
% Standard Gaussian Mixture Model fit to the data. I got this Matlab code from Todd Leen.
% 
% pnormal.m
% used by mixturefit.m 
%
% plotgmmpdf.m
% A (very stupid) function that plots the GMM pdf. Only works in 2 dimensions. 
% It inputs the GMM parameters along with the boundaries and the sampling density of the pdf to be plotted.
% 
% pc_projection_gmm.m
% This is the main routine that projects any point x into the feature space 
% onto the principal surface of target dimensionality. 
% - Inputs the GMM parameters, and x. 
% - Outputs x* - principal curve projection of x.
%
% pc_projection_multidim.m
% This is the main pc projection routine based on KDE. inputs the data
% samples and the kernel size, outputs the principal curve projection of
% the data
%
% gmmpdf.m
% - Inputs the GMM parameters and x
% - Outputs p(x), g(x), H(x), SI(x) (and also optionally another vector for an
% intermediate step that I use in the projection)
% 
%%%%%%%%%%%%%%%

clear all
close all

% generating data
% dataflag = 1   A single Gaussian in 2D
% dataflag = 2   Gaussian mixture with 3 components in 2D (*-type intersection)
% dataflag = 3   Gaussian mixture with 3 components in 2D (T-type intersection)

n=501;
data = zeros(3,n);
r1 = randn(1,n);
r2 = randn(1,n);
r3 = randn(1,n);
t = linspace(0,10*pi,n);
x = 5*sin(t);
y = 5*cos(t);
z = t;
data(1,:) = x+0.1*r1;
data(2,:) = y+0.1*r2;
data(3,:) = z+0.5*r3;

plot3(data(1,:),data(2,:),data(3,:),'p'),axis equal;grid;

% fitting GMM to the data 

% assuming we know the number of components
%if dataflag==1,ncomp=1;else,ncomp=3;end

%constraint=0; %which means we are looking for full covariance
%eps=1e-6; %default value
%[alphas,means,mcovs]=mixturefit(data',ncomp,constraint,eps);

% let's see how it looks like
%scale=trace(cov(data'));
%xmin=min(data(1,:))-scale/10;ymin=min(data(2,:))-scale/10;
%xmax=max(data(1,:))+scale/10;ymax=max(data(2,:))+scale/10;
%density=10;
%plotgmmpdf(xmin,xmax,ymin,ymax,density,alphas,means,mcovs);

% now project the data points onto the principal curve
targetdim = 1; % the target dimensionality is 1 for principal curve. 
               % for any principal surface of different dimensionality just 
               % change this
data_in = data;% initialize on a grid, or onto the data points themselves
%pc_projection = pc_projection_GMM(means,alphas,mcovs,targetdim,data_in)

%figure,plot(data(1,:),data(2,:),'.');hold on,
%title('PC projection with GMM pdf')
%plot(pc_projection(1,:),pc_projection(2,:),'.r');axis equal;

% now let's see what happens with a KDE... 
% targetdim = 1 as before, we don't want to find a principal surface etc.
% but the principal curve.
% 
kernel_sigma = 1.5; % play with this parameter... I handpicked the one over
                    % here but there are many ways to select a good kernel 
                    % bandwidth. See my thesis for references...

pc_projection=pc_project_multidim(data,data_in,kernel_sigma,targetdim)
pc_projection = pc_projection';
figure,plot3(data(1,:),data(2,:),data(3,:),'.');hold on,
title('PC projection with KDE pdf')
plot3(pc_projection(1,:),pc_projection(2,:),pc_projection(3,:),'.r');axis equal;grid;

keyboard;