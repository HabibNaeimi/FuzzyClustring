clc
clear all
close all
%% OFFFLINE & ONLINE CLUSTRING ALGORITHEM.
%   This programm developed as a course project for Fuzzu systems course by
%    Habibollah Naeimi.
disp(' OFFFLINE & ONLINE CLUSTRING ALGORITHEM.');
disp('  This programm develops as a course project for Fuzzu systems course')
disp('   by Habibollah Naeimi.')
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')

%% 1st Part: Parameter Setting.
%%  Parameters Initiating.

Sigma = 2;                           % Gaussian Membership function spread.
%Sigma = input(' Please Enter the spread of gaussian MF:');
Radius = 0.2;                        % Radius of clusters.
%Radius = input(' Please Enter the Radius of Clusters:');

DataPairNu = 100;                    % Number of Data Pairs.
%DataPairNu = input(' Please Enter the Number of Data Pairs:');
SampleNum = 500;                     % Number of Samples.
%SampleNum = input(' Please Enter the Number of samples:');

InpNum = 2;                          % Number of Inputs.

Clustring_Mode = 1;                  
disp(' Clustring Modes: ');     
disp(' Offline: First System Identifying then Controllimg. Enter 1.');
disp(' Online: System Identifying and Controlling simultaneously. Enter 2.');
disp(' ');
Clustring_Mode = input(' Please Enter The Clustring Mode:');
disp(' ');
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');
%% Offline Mode.
if Clustring_Mode==1
    tic                              % Starting Timer.
%% 1st Part: Data Generation.

    y = zeros(1,DataPairNu+InpNum);
    Pairs = zeros(DataPairNu,InpNum+1);
    y(1:2) = [0.1 0.2];              % Initialing first two of y.
        
    for k=3:DataPairNu+InpNum        % Generating Data.
        
        r = sin(2*pi*k/25);            
        g(k) = y(k-1)*y(k-2)*(y(k-1)+2.5)/(1+y(k-1)^2+y(k-2)^2);
        y(k) = r+g(k);
            
    end        
    for i=1:DataPairNu
        Pairs(i,:)=y(i:i+InpNum);    % Generating Data Pairs.
    end
        
%% 2nd Part: Rule Generation.

    n = 1;
    DataPairNu = size(Pairs,1);
    x_centr(1,:) = Pairs(1,1:end-1); % Eatablishing First Cluster Center.
    A(1) = Pairs(1,end);             % and its A and B.
    B(1) = 1;
    
    for p=2:DataPairNu
    
        x_x = repmat(Pairs(p,1:end-1),size(x_centr,1),1);
        FinalVAL = abs(x_centr-repmat(Pairs(p,1:end-1),size(x_centr,1),1));
        DISTNS = max(FinalVAL,[],2);
        Indx = find(DISTNS<=Radius);
        
        if isempty(Indx)
           x_centr = [x_centr;Pairs(p,1:end-1)];
           A = [A;Pairs(p,end)];
           B = [B;1];
           n = n+1;
        else
           
           A(Indx(1),:) = A(Indx(1),:)+Pairs(p,end);   % Updating Clusters.
           B(Indx(1),:) = B(Indx(1),:)+1;
        end
    end
    
    disp(' Number of Clusters:');
    disp(n);                         % Number of Clusters.        
    
%% 3rd Part: Plotting Clusters.

    figure;                          % Plotting Datas.
    plot(Pairs(:,end));
    title('Just Training Datas.');
    
    figure;                          % Plotting Clusters in 2 dim.
    plot(x_centr(:,1),x_centr(:,2),'*r');
    title('Clusters in two Dim.');
    
    figure;                          % Plotting Clusters in 3 dim.
    Y = A./B;
    plot3(x_centr(:,1),x_centr(:,2),Y,'*r');
    grid on
    title('Clusters in three Dim.');
        
%% 4th Part: Result.        
        
    y_m = zeros(1,SampleNum);        % Initialing.
    f = zeros(1,SampleNum);
    g = zeros(1,SampleNum);
    y_Es = zeros(1,SampleNum);
    y_Es(1:3) = [0.1 0.2 0.3];
        
    for k=3:SampleNum
        f(k) = TheF(x_centr,A,B,y_Es(k-InpNum:k),Sigma);
        g(k) = y_Es(k-1)*y_Es(k-2)*(y_Es(k-1)+2.5)/(1+y_Es(k-1)^2+y_Es(k-2)^2);

        u = -f(k)+0.6*y_Es(k-1)+0.2*y_Es(k-2)+sin(2*pi*(k-1)/25);
        y_Es(k) = 0.6*y_Es(k-1)+0.2*y_Es(k-2)+sin(2*pi*(k-1)/25);
        
        y_m(k) = u+g(k);                      
    end  
        
%% 5th Part: Error and Final Plotting.

    MSE=mse(y_Es-y_m)               % Mean Squared Error.
    
    figure;                         % Plotting Results.
    plot(y_Es)         
    hold on
    plot(y_m,'r')
    legend('Estimated Result','Current or Real Value');
    toc                             % Stopping Timer.
        
%% Online Mode.
elseif Clustring_Mode==2;
       tic                          % Start Timer.
%% 1st Part: Data and Rule Generation.

    y(1:2)=[0.5 1];
    R = sin(2*pi*[1:DataPairNu]/25);
    x_centr(:,1) = [y(1);y(2)];     % Establishing First Cluster Center.
    A(1) = (y(1)*y(2)*(y(2)+2.5))/(1+y(1)^2+y(2)^2);
    B(1) = 1;
    M = 1;
    
    for k = 2:DataPairNu-1
        f = OnlineF(x_centr,y(k-1),y(k),Sigma,A,B);
        p = (y(k-1)*y(k)*(y(k)+2.5))/(1+y(k-1)^2+y(k)^2);
        y(k+1) = p-f+0.6*y(k)+0.2*y(k-1)+R(k);
    
        for i= 1:M
            Dist(i) = (sqrt((y(k-1)-x_centr(1,i))^2+(y(k)-x_centr(2,i))^2)-Radius);
        end
    
        [disValue,index] = min(Dist);
        
        if disValue < 0
                                    % Updating Cluster.
            A(index) = A(index)+(y(k-1)*y(k)*(y(k)+2.5))/(1+y(k-1)^2+y(k)^2);
            B(index) = B(index)+1;
        else
            
            M = M+1;                % Creating New Clusters.
            x_centr(:,M) = [y(k-1);y(k)];
            A(M) = (y(k-1)*y(k)*(y(k)+2.5))/(1+y(k-1)^2+y(k)^2);
            B(M) = 1;
        end
    end
%%  2nd Part: Result Calculating.

    ym(1:2) = y(1:2);
    for t=2:DataPairNu-1
        ym(t+1) = 0.6*ym(t)+0.2*ym(t-1)+R(t);
    end
    
%%  3rd Part: Error and Plotting.

    Error = y(3:DataPairNu)-ym(3:DataPairNu);
    MSE = mse(Error)                % Mean Squared Error.
    
    disp(' Number of Clusters:');
    disp(M);
    
    figure                          % Plotting Result.
    plot(1:DataPairNu,ym,1:DataPairNu,y,'r')
    title('Online Clustring.')
    legend('Y',' Estimated Y')
    
    figure(2)                       % Plotting Clusters.
    scatter(x_centr(1,:),x_centr(2,:),'m','filled')
    title('The Clusters.')
    toc                             % Stopping Timer.

    
else
    disp(' Wrong Answer! Start Again.');
end

