close all;
clearvars;
clc;

projectRoot=fileparts(mfilename('fullpath'));
addpath(projectRoot);

%% General configuration
params.fonction='rosenbrock';
params.b=2;
Xoptimal=[1;1];

options.tolX=1e-8;
options.tolF=1e-8;
options.tolG=1e-8;
options.maxiter=1000;
options.maxBacktracking=200;
options.divergenceFactor=1e12;

% Armijo backtracking parameters from the assignment
options.initialStep=1;
options.armijoC=1e-4;
options.backtrackFactor=0.75;

% Fixed-step and Levenberg-Marquardt parameters
options.pas='fixe';
options.step=1;
options.lmMuInitial=1e-3;
options.lmMuMin=1e-12;
options.lmMuMax=1e12;

%% 1.1 - Numerical verification of the analytical derivatives
derivativeChecks=verify_derivatives(params);
disp('Finite-difference derivative checks:');
disp(derivativeChecks);

%% 1.2 - Rosenbrock objective visualization
x1=-4:0.05:4;
x2=-5:0.05:20;
[X1,X2]=meshgrid(x1,x2);
Fgrille=zeros(size(X1));

for k=1:numel(X1)
    Fgrille(k)=rosenbrock([X1(k);X2(k)],params);
end

figure('Name','Rosenbrock surface','Color','w');
meshc(X1,X2,Fgrille);
xlabel('x_1');
ylabel('x_2');
zlabel('f(x)');
title('Rosenbrock objective, b=2');
grid on;

figure('Name','Rosenbrock contours','Color','w');
contour(X1,X2,Fgrille,100);
hold on;
plot(Xoptimal(1),Xoptimal(2),'rp', ...
    'MarkerFaceColor','r','MarkerSize',12);
xlabel('x_1');
ylabel('x_2');
title('Rosenbrock objective contours');
legend('Level sets','Analytical minimizer (1,1)','Location','best');
grid on;

%% 1.3.1 - Fixed-step gradient descent
X0=[0;0];
pasFixes=[1,1e-2,1e-4];
resultatsFixes=cell(length(pasFixes),1);
lignesFixes=cell(length(pasFixes),7);

options.method='gradient';
options.pas='fixe';

for k=1:length(pasFixes)
    options.step=pasFixes(k);
    [xh,result,xval]=optimdescent(params,options,X0);

    resultatsFixes{k}=struct('xh',xh,'result',result,'xval',xval);
    lignesFixes(k,:)={pasFixes(k),xh(1),xh(2),result.crit(end), ...
        result.grad(end),result.iter,result.stop};
end

tableFixe=cell2table(lignesFixes,'VariableNames', ...
    {'Step','x1','x2','Objective','GradientNorm','Iterations','Stop'});
disp('Fixed-step gradient descent:');
disp(tableFixe);

figure('Name','Fixed-step convergence','Color','w');
tiledlayout(1,length(pasFixes),'TileSpacing','compact');

for k=1:length(pasFixes)
    nexttile;
    valeurs=resultatsFixes{k}.result.crit;
    semilogy(0:length(valeurs)-1,max(valeurs,realmin),'LineWidth',1.5);
    xlabel('Iteration');
    ylabel('f(x)');
    title(sprintf('Step = %g (%s)',pasFixes(k), ...
        resultatsFixes{k}.result.stop));
    grid on;
end
sgtitle('Sensitivity to the fixed step');

figure('Name','Fixed-step trajectories','Color','w');
tiledlayout(1,length(pasFixes),'TileSpacing','compact');

for k=1:length(pasFixes)
    nexttile;
    contour(X1,X2,Fgrille,50);
    hold on;

    valeursX=resultatsFixes{k}.xval;
    valides=all(isfinite(valeursX),1) ...
        & valeursX(1,:) >= -4 & valeursX(1,:) <= 4 ...
        & valeursX(2,:) >= -5 & valeursX(2,:) <= 20;
    trajectoire=valeursX;
    trajectoire(:,~valides)=NaN;

    plot(trajectoire(1,:),trajectoire(2,:),'-o', ...
        'LineWidth',1.2,'MarkerSize',3);
    plot(1,1,'rp','MarkerFaceColor','r','MarkerSize',10);
    xlim([-4,4]);
    ylim([-5,20]);
    xlabel('x_1');
    ylabel('x_2');
    title(sprintf('Step = %g',pasFixes(k)));
    grid on;
end

%% 1.3.2 - Gradient descent with Armijo backtracking
options.method='gradient';
options.pas='variable';

[xhArmijo,resultArmijo,xvalArmijo]=optimdescent(params,options,X0);

fprintf(['\nGradient with Armijo: x=(%.8g, %.8g), f=%.3e, ', ...
    '|g|=%.3e, iter=%d, %s\n'],xhArmijo(1),xhArmijo(2), ...
    resultArmijo.crit(end),resultArmijo.grad(end), ...
    resultArmijo.iter,resultArmijo.stop);

figure('Name','Armijo gradient descent','Color','w');
tiledlayout(1,3,'TileSpacing','compact');

nexttile;
semilogy(0:length(resultArmijo.crit)-1, ...
    max(resultArmijo.crit,realmin),'LineWidth',1.5);
xlabel('Iteration');
ylabel('f(x)');
title('Objective');
grid on;

nexttile;
semilogy(0:length(resultArmijo.grad)-1, ...
    max(resultArmijo.grad,realmin),'LineWidth',1.5);
xlabel('Iteration');
ylabel('||grad f(x)||_2');
title('Gradient norm');
grid on;

nexttile;
contour(X1,X2,Fgrille,50);
hold on;
plot(xvalArmijo(1,:),xvalArmijo(2,:),'-o','MarkerSize',3);
plot(1,1,'rp','MarkerFaceColor','r','MarkerSize',10);
xlim([-4,4]);
ylim([-5,20]);
xlabel('x_1');
ylabel('x_2');
title('Trajectory');
grid on;

%% 1.3.3-1.3.6 - Comparison of the seven configurations
configurations=struct( ...
    'nom',{'Gradient Armijo','Fletcher-Reeves','Newton unit step', ...
           'Newton Armijo','BFGS','Gauss-Newton','Levenberg-Marquardt'}, ...
    'method',{'gradient','gradient conjuge','newton','newton', ...
              'quasi-newton','gauss-newton','levenberg-marquardt'}, ...
    'pas',{'variable','variable','fixe','variable', ...
           'variable','variable','variable'}, ...
    'step',{1,1,1,1,1,1,1});

pointsInitiaux={[-2;10],[4;-5]};
comparaisons=cell(length(pointsInitiaux),length(configurations));
lignesComparaison={};
couleurs=lines(length(configurations));

for ip=1:length(pointsInitiaux)
    for im=1:length(configurations)
        options.method=configurations(im).method;
        options.pas=configurations(im).pas;
        options.step=configurations(im).step;

        [xh,result,xval]=optimdescent(params,options,pointsInitiaux{ip});
        comparaisons{ip,im}=struct('xh',xh,'result',result,'xval',xval);

        lignesComparaison(end+1,:)={ ...
            sprintf('[%g;%g]',pointsInitiaux{ip}(1),pointsInitiaux{ip}(2)), ...
            configurations(im).nom,xh(1),xh(2),result.crit(end), ...
            result.grad(end),result.iter,result.temp,result.stop}; %#ok<SAGROW>
    end

    figure('Name',sprintf('Method comparison - start %d',ip),'Color','w');
    tiledlayout(1,2,'TileSpacing','compact');

    nexttile;
    for im=1:length(configurations)
        valeurs=comparaisons{ip,im}.result.crit;
        semilogy(0:length(valeurs)-1,max(valeurs,realmin), ...
            'LineWidth',1.3,'Color',couleurs(im,:));
        hold on;
    end
    xlabel('Iteration');
    ylabel('f(x)');
    title(sprintf('Objective, x_0=[%g;%g]', ...
        pointsInitiaux{ip}(1),pointsInitiaux{ip}(2)));
    legend({configurations.nom},'Location','best');
    grid on;

    nexttile;
    for im=1:length(configurations)
        valeurs=comparaisons{ip,im}.result.grad;
        semilogy(0:length(valeurs)-1,max(valeurs,realmin), ...
            'LineWidth',1.3,'Color',couleurs(im,:));
        hold on;
    end
    xlabel('Iteration');
    ylabel('||grad f(x)||_2');
    title(sprintf('Gradient norm, x_0=[%g;%g]', ...
        pointsInitiaux{ip}(1),pointsInitiaux{ip}(2)));
    legend({configurations.nom},'Location','best');
    grid on;
end

tableComparaison=cell2table(lignesComparaison,'VariableNames', ...
    {'InitialPoint','Method','x1','x2','Objective','GradientNorm', ...
     'Iterations','TimeSeconds','Stop'});

disp('Comparison of optimization methods:');
disp(tableComparaison);
