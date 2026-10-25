function checks=verify_derivatives(params)
%VERIFY_DERIVATIVES Controle g, H et J par differences finies centrees.

    dimensions=[2,3,5];
    nbTests=2*length(dimensions);
    lignes=cell(nbTests,5);
    ligne=0;

    for n=dimensions
        points={linspace(-1.2,1.1,n)',sin((1:n)')};

        for k=1:length(points)
            x=points{k};
            [~,g,H,r,J]=rosenbrock(x,params);

            pasG=1e-6;
            pasH=1e-5;
            gFD=zeros(n,1);
            HFD=zeros(n,n);
            JFD=zeros(length(r),n);

            for j=1:n
                e=zeros(n,1);
                e(j)=1;

                [fPlus,gPlus,~,rPlus,~]= ...
                    rosenbrock(x+pasG*e,params);
                [fMoins,gMoins,~,rMoins,~]= ...
                    rosenbrock(x-pasG*e,params);

                gFD(j)=(fPlus-fMoins)/(2*pasG);
                JFD(:,j)=(rPlus-rMoins)/(2*pasG);

                [~,gPlusH]=rosenbrock(x+pasH*e,params);
                [~,gMoinsH]=rosenbrock(x-pasH*e,params);
                HFD(:,j)=(gPlusH-gMoinsH)/(2*pasH);
            end

            erreurG=norm(g-gFD)/max(1,norm(g));
            erreurH=norm(H-HFD,'fro')/max(1,norm(H,'fro'));
            erreurJ=norm(J-JFD,'fro')/max(1,norm(J,'fro'));
            erreurMC=norm(g-2*J'*r)/max(1,norm(g));

            ligne=ligne+1;
            lignes(ligne,:)={n,k,erreurG,erreurH,max(erreurJ,erreurMC)};
        end
    end

    checks=cell2table(lignes,'VariableNames', ...
        {'Dimension','Point','GradientError','HessianError','LeastSquaresError'});

    erreurs=[checks.GradientError;checks.HessianError;checks.LeastSquaresError];
    if any(~isfinite(erreurs)) || max(erreurs) > 1e-5
        error('verify_derivatives:Failure', ...
            'La verification numerique des derivees a echoue.');
    end
end
