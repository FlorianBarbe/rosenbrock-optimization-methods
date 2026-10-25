%EXPORT_FIGURES Regenerates and exports the figures used by the report.

main;

figureDir=fullfile(projectRoot,'figures');
if ~exist(figureDir,'dir')
    mkdir(figureDir);
end

exports={ ...
    'Rosenbrock surface','partie1_01_rosenbrock_surface.png'; ...
    'Rosenbrock contours','partie1_02_rosenbrock_contours.png'; ...
    'Fixed-step convergence','partie1_03_rosenbrock_pas_fixes.png'; ...
    'Fixed-step trajectories','partie1_04_rosenbrock_trajectoires_pas_fixes.png'; ...
    'Armijo gradient descent','partie1_05_rosenbrock_armijo.png'; ...
    'Method comparison - start 1','partie1_06_comparaison_rosenbrock_point_1.png'; ...
    'Method comparison - start 2','partie1_07_comparaison_rosenbrock_point_2.png'};

for k=1:size(exports,1)
    figureHandle=findall(groot,'Type','figure','Name',exports{k,1});
    if isempty(figureHandle)
        error('export_figures:MissingFigure', ...
            'Figure not found: %s',exports{k,1});
    end

    set(figureHandle(1),'Color','w');
    exportgraphics(figureHandle(1),fullfile(figureDir,exports{k,2}), ...
        'Resolution',220);
end

writetable(tableFixe,fullfile(figureDir,'rosenbrock_pas_fixes.csv'));
writetable(tableComparaison,fullfile(figureDir,'rosenbrock_comparaison.csv'));
writetable(derivativeChecks,fullfile(figureDir,'rosenbrock_derivative_checks.csv'));

fprintf('Figures and tables exported to %s\n',figureDir);
