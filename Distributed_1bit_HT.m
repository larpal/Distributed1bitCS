%% Simulation of distributed one-bit compressive sensing
% This code runs the simulations for the paper "Analysis of
% Hard-Thresholding for Distributed Compressed Sensing with One-Bit
% Measurements" by J. Maly and L. Palzer.
%%%%%%%
close all
clear

       dim = 100;           % Signal dimension
sparsity   = .025:.025:.5;  % Fraction of nonzero coefficients - we used .005:.005:0.5 in the paper
rate       = .5:.5:3;       % Nr. of measurements/dim - we used .005:.005:3 in the paper 
nUser      = [2,5,20];      % Number of users 
nsim       = 50;            % Number of simulations for every parameter - we used 500 in the paper

 Err_joint = zeros(nsim,length(rate),length(sparsity),length(nUser));   % array to store errors for distributed CS
  Err_sing = zeros(nsim,length(rate),length(sparsity));                 % array to store errors for single-user CS
%
fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,nsim) '\n\n']);

parfor iter = 1:nsim
        err_d1 = ones(length(rate),length(sparsity),length(nUser));    % tmp error string within parfor
         err_s = ones(length(rate),length(sparsity));                  % tmp error string within parfor
    for p = 1:length(sparsity)  
        for r = 1:length(rate)
%% Single user 1-bit CS by hard thresholding  
            b = [ones(round(dim*sparsity(p)),1) ; zeros(dim-round(dim*sparsity(p)),1)];
            b = b(randperm(dim));                                               % create sparsity patern
            G = randn(dim,1);                                                   % create Gaussian nonzero coefficients
            X = b.*G;                                                           % create signal
            X = X./norms(X,2,1);                                                % scale signal
         meas = ceil(rate(r)*dim);
                   A = randn(meas,dim)*sqrt(pi/(2*meas^2));                     % scale matrix
                   Y = sign(A * X);                                             % 1-bit measurements
                  z1 = HT_1Bit( A(1:meas,:),Y(:,1),round(dim*sparsity(p)) );    % recover signal
          err_s(r,p) = norm(X(:,1)-z1);                                         % compute error

%% Distributed 1-bit CS via joint hard thresholding          
            for u = 1:length(nUser)
                b = [ones(round(dim*sparsity(p)),1) ; zeros(dim-round(dim*sparsity(p)),1)];
                b = b(randperm(dim));                                           % create sparsity patern
                G = randn(dim,nUser(u));                                        % create Gaussian nonzero coefficients
                X = b.*G;                                                       % create signal
                X = X./norms(X,2,1)/sqrt(nUser(u));                             % scale signal
                meas = ceil(rate(r)*dim);
                   A = randn(meas*nUser(u),dim)*sqrt(pi/(2*nUser(u)*meas^2));   % scale matrix
                   Y = [];
                   for k = 1:nUser(u)
                       Y = [Y, sign(A( (k-1)*meas+1:k*meas,: ) * X(:,k))];      % 1-bit measurements
                   end
                   Z1 = joint_HT_1Bit_distributed( A,Y,dim,meas,nUser(u),round(dim*sparsity(p)) );  %recover signal
        err_d1(r,p,u) = norm(X-Z1,'fro');                                       % compute error
            end
        end
    end
    Err_joint(iter,:,:,:) = err_d1;
    Err_sing(iter,:,:) = err_s;
    fprintf('\b|\n');
end
%
%% Save
 Filename = sprintf('HT_1Bit_%s.mat', datestr(now,'mm-dd-yyyy HH-MM'));
 save(Filename,'rate','dim',...
     'sparsity','nsim','Err_joint','Err_sing','nUser');

%% Loglog Plot
   lw = 2;
figure
loglog(rate,median(Err_sing(:,:,2)),'LineWidth',lw)
hold on
loglog(rate,median(Err_joint(:,:,2,1)),'LineWidth',lw)
loglog(rate,median(Err_joint(:,:,2,2)),'LineWidth',lw)
loglog(rate,median(Err_joint(:,:,2,3)),'LineWidth',lw)
grid on
legend('L = 1','L = 2','L = 5','L = 20')
xlabel('measurement rate m/n')
ylabel('average X-Xhat_F')

%% Phase Plot
D_av_joint  = squeeze(mean(Err_joint,1));
 D_av_sing  = squeeze(mean(Err_sing,1));
 d_recovery = 2/3;        % recovery threshold
%
figure
hold on
subplot(2,3,1); 
imagesc(sparsity,rate,(squeeze(D_av_sing)));
set(gca,'YDir','normal')
titlestr = sprintf('L=%d user',1);
title(titlestr)
ylabel('measurement rate m/n')
xlabel('sparsity s/n')
hold on
contour(sparsity,rate,(squeeze(D_av_sing)),[d_recovery d_recovery],'r')
caxis([0 1])
%
for user = 1:length(nUser)
subplot(2,3,1+user+(user>1)); 
imagesc(sparsity,rate,(squeeze(D_av_joint(:,:,user))));
set(gca,'YDir','normal')
titlestr = sprintf('L=%d users',nUser(user));
title(titlestr)
ylabel('measurement rate m/n')
xlabel('sparsity s/n')
hold on
contour(sparsity,rate,(squeeze(D_av_joint(:,:,user))),[d_recovery d_recovery],'r')
caxis([0 1])
end
subplot(2,3,[3 6]); caxis([0 1]);colorbar('West');set(gca,'vis','off');


