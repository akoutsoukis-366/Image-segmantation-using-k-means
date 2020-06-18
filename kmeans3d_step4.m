% --- initialize Matlab
clc;
disp('-----');
clear;
close all ;

% --- initialize the random generator
% --- for debugging purposes (creates always the same "random" numbers)
%rng('default');

% --- select a file (image) to load the 3D data points (pixels)
Im=imread('butterfly.jpg');

% --- show the original image
figure;
image(Im);
axis image;

% --- define k
K=5;

% --- get the dimensions of the image and the number of color bands
[R,Co,Bands]=size(Im);

% ---calc the overall number of pixels
N=R*Co;

% --- reshape the image data from RxCx3 (Rows x Columns x ColorBands) into an Nx3 array (r,g,b triplets)
% --- where N=RxC, that holds the (r,g,b) color values for all the N pixels
X=double(reshape(Im,N,Bands));

C=[];
while size(unique(C,'rows'),1)~=K% ensure that there are k DIFFERENT centers
    % --- find k random unique pixels (i.e. with different rgb values)
    Idx=round(rand(K,1)*N);    
    % --- set these pixels as initial cluster centers
    C=X(Idx,:);
end
old_C=zeros(K,3);
new_C=ones(K,3);

Iter=0;
while Iter~=50;
Iter=Iter+1;
    % --- find the distance of all N points from the K clusters.
    % --- the resulting array Dist has dimensions [NxK].
    % --- for example Dist(7,2) is the distance of the 7th data point from the 2nd center.
    for i=1:K
        Dist(:,i)=sqrt(sum((X-repmat(C(i,:),N,1)).^2,2));
       old_C=C; 
      
    end
     [MinValue,T]=min(Dist,[],2);
    % --- for each one of the N data points, find the minimum of the K distances (MinValue)
    % --- and the corresponding cluster (T)
    %[MinValue,T]=min(Dist,[],2);
    
    % --- update each cluster's center by the points that are assigned to it
    for i=1:K
        if sum(T==i)>0
            C(i,:)=(mean(X(T==i,:),1));
           new_C=C;
        end
    end

   % ?? ????? ???? ???????????? ??? ??? ?? ? ????????
m=mean(X,1);

for i=1:K
     if sum(T==i)>0
   E_distance=sqrt(sum(C(i,:)-m));
   E_distance_1=E_distance.^2;
   s=size((X(T==i)),1);
   b=s*E_distance;
   B=b/(1-K);
     end;
   
end
W=sum(MinValue);

W=W/(size(N,1)-K);

end



cla;
T1=T(1:200:end,1);
% --- show the pixels (NOT ALL, just every 500-th pixel) as points in the RGB space
scatter3(X(1:200:end,1),X(1:200:end,2),X(1:200:end,3),20,T1,'filled');
hold on;


% --- set drawing parameters
axis equal;
box on;
view([-17 30]);
xlabel('Red','color','r');
ylabel('Green','color','g');
zlabel('Blue','color','b');
colormap(C/255);


% --- show the cluster centers with k different colors
scatter3(C(:,1),C(:,2),C(:,3),100,(1:K)','o','filled','MarkerEdgeColor','k');
hold on;
scatter3(m(:,1),m(:,2),m(:,3),200,'o','filled','MarkerEdgeColor','k');
 title(['Iteration ' num2str(Iter) ],'Color','r');
 drawnow;



Im2=(T(:,1));

Im2=reshape(Im2,[R,Co]);
figure;
image(Im2);
colormap(C/255);
axis image;

