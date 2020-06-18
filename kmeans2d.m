% --- initialize
disp('-----');
clc;
clear;
close all;

% --- select a file to load the 2D data points
A=load('points2d.txt');

% --- select a file to load the 2D initial clulster centers
C=load('centers2d.txt');

% --- get the number of data points
N=size(A,1);

% --- get the number of clusters K
K=size(C,1);

% --- define the cluster colors in the figure according to number K
Colors=hsv(max(3,K))
colormap(Colors);

% --- plot the data points and the original centers
hold on;
scatter(A(:,1),A(:,2),10,'k','filled');
scatter(C(:,1),C(:,2),50,(1:K)','filled');


% --- set drawing parameters
set(gca,'ticklength',[0 0]);
box on;
axis image;
axis([0 max([A(:,1);C(:,1)])+20 0 max([A(:,2);C(:,2)])+20]);

title('Initialization - Press a key to start...','FontSize',14,'Color','r');
pause;

% --- make the first 5 Iterations
for Iter=1:5       

    % --- find the distance of all N points from the K clusters.
    % --- the resulting array Dist has dimensions [NxK].
    % --- for example Dist(7,2) is the distance of the 7th data point from the 2nd center.
    for i=1:K
        Dist(:,i)=sqrt(sum((A-repmat(C(i,:),N,1)).^2,2));
    end
    
    % --- for each one of the N data points, find the minimum of the K distances (MinValue)
    % --- and the corresponding cluster (T)
    [MinValue,T]=min(Dist,[],2);
    
    % --- update each cluster's center by the points that are assigned to it
    for i=1:K
        if sum(T==i)>0
            C(i,:)=(mean(A(T==i,:),1));
        end
    end
    
    % --- color the data points according to the cluster they belong
    % --- and show the updated centers position
    cla;
    scatter(A(:,1),A(:,2),10,T,'filled');
    scatter(C(:,1),C(:,2),50,(1:K)','filled');    
    title(['Iteration ' num2str(Iter) '/5'],'FontSize',14,'Color','r');
    pause(1);
end


