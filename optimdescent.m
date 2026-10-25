function [xh,result,xval]=optimdescent(params,options,x0)
%OPTIMDESCENT Methodes de descente avec recherche lineaire.
%   [xh,result,xval]=optimdescent(params,options,x0) minimise la fonction
%   de Rosenbrock en utilisant la methode choisie dans options.method.

    tempsDebut=tic;

    X=x0(:);
    [F,G,H,r,J]=rosenbrock(X,params);

    xh=X;
    xval=X;

    result.iter=0;
    result.crit=F;
    result.grad=norm(G);
    result.step=[];
    result.stop='MaxIter';
    result.mu=[];

    n=length(X);
    critereInitial=max(1,abs(F));

    if norm(G) <= options.tolG
        result.stop='TolG';
        result.temp=toc(tempsDebut);
        return
    end

    switch options.method
        case 'gradient'
            switch options.pas

%%============ DESCENTE DE GRADIENT A PAS FIXE ============
                case 'fixe'
                    alpha=options.step;

                    for i=1:options.maxiter
                        direction=-G;
                        Xnew=X+alpha*direction;
                        [Fnew,Gnew,~,~,~]=rosenbrock(Xnew,params);

                        if estDivergent(Fnew,Gnew,critereInitial,options)
                            result.stop='Divergence';
                            result.iter=i-1;
                            break
                        end

                        [continuer,result]=testsortie( ...
                            X,Xnew,F,Fnew,Gnew,i,result,options);

                        xval(:,i+1)=Xnew;
                        result.crit(i+1)=Fnew;
                        result.grad(i+1)=norm(Gnew);
                        result.step(i)=alpha;
                        result.iter=i;

                        X=Xnew;
                        F=Fnew;
                        G=Gnew;

                        if ~continuer
                            break
                        end
                    end

%%============ DESCENTE DE GRADIENT AVEC ARMIJO ============
                case 'variable'
                    for i=1:options.maxiter
                        direction=-G;

                        [pasAccepte,alpha,Xnew,Fnew,Gnew,~,~,~]= ...
                            rechercheArmijo(X,F,G,direction,params,options);

                        if ~pasAccepte
                            result.stop='RecherchePas';
                            result.iter=i-1;
                            break
                        end

                        [continuer,result]=testsortie( ...
                            X,Xnew,F,Fnew,Gnew,i,result,options);

                        xval(:,i+1)=Xnew;
                        result.crit(i+1)=Fnew;
                        result.grad(i+1)=norm(Gnew);
                        result.step(i)=alpha;
                        result.iter=i;

                        X=Xnew;
                        F=Fnew;
                        G=Gnew;

                        if ~continuer
                            break
                        end
                    end

                otherwise
                    error('optimdescent:Pas','Type de pas inconnu: %s',options.pas);
            end

%%============ GRADIENT CONJUGUE NON LINEAIRE ============
        case 'gradient conjuge'
            direction=-G;

            for i=1:options.maxiter
                if G'*direction >= 0 || any(~isfinite(direction))
                    direction=-G;
                end

                [pasAccepte,alpha,Xnew,Fnew,Gnew,~,~,~]= ...
                    rechercheArmijo(X,F,G,direction,params,options);

                if ~pasAccepte
                    result.stop='RecherchePas';
                    result.iter=i-1;
                    break
                end

                [continuer,result]=testsortie( ...
                    X,Xnew,F,Fnew,Gnew,i,result,options);

                xval(:,i+1)=Xnew;
                result.crit(i+1)=Fnew;
                result.grad(i+1)=norm(Gnew);
                result.step(i)=alpha;
                result.iter=i;

                if continuer
                    betaFR=(Gnew'*Gnew)/max(G'*G,eps);
                    directionNew=-Gnew+betaFR*direction;

                    if Gnew'*directionNew >= 0 || any(~isfinite(directionNew))
                        directionNew=-Gnew;
                    end
                end

                X=Xnew;
                F=Fnew;
                G=Gnew;

                if ~continuer
                    break
                end

                direction=directionNew;
            end

%%============ METHODE DE NEWTON ============
        case 'newton'
            for i=1:options.maxiter
                direction=directionNewton(H,G);

                if G'*direction >= 0 || any(~isfinite(direction))
                    direction=-G;
                end

                if strcmp(options.pas,'fixe')
                    alpha=options.step;
                    Xnew=X+alpha*direction;
                    [Fnew,Gnew,Hnew,~,~]=rosenbrock(Xnew,params);
                    pasAccepte=~estDivergent( ...
                        Fnew,Gnew,critereInitial,options);

                    if ~pasAccepte
                        result.stop='Divergence';
                        result.iter=i-1;
                        break
                    end
                else
                    [pasAccepte,alpha,Xnew,Fnew,Gnew,Hnew,~,~]= ...
                        rechercheArmijo(X,F,G,direction,params,options);

                    if ~pasAccepte
                        result.stop='RecherchePas';
                        result.iter=i-1;
                        break
                    end
                end

                [continuer,result]=testsortie( ...
                    X,Xnew,F,Fnew,Gnew,i,result,options);

                xval(:,i+1)=Xnew;
                result.crit(i+1)=Fnew;
                result.grad(i+1)=norm(Gnew);
                result.step(i)=alpha;
                result.iter=i;

                X=Xnew;
                F=Fnew;
                G=Gnew;
                H=Hnew;

                if ~continuer
                    break
                end
            end

%%============ METHODE QUASI-NEWTON BFGS ============
        case 'quasi-newton'
            I=eye(n);
            B=I;

            for i=1:options.maxiter
                direction=-B*G;

                if G'*direction >= 0 || any(~isfinite(direction))
                    B=I;
                    direction=-G;
                end

                [pasAccepte,alpha,Xnew,Fnew,Gnew,~,~,~]= ...
                    rechercheArmijo(X,F,G,direction,params,options);

                if ~pasAccepte
                    result.stop='RecherchePas';
                    result.iter=i-1;
                    break
                end

                [continuer,result]=testsortie( ...
                    X,Xnew,F,Fnew,Gnew,i,result,options);

                if continuer
                    s=Xnew-X;
                    y=Gnew-G;
                    courbure=y'*s;
                    seuil=1e-12*norm(s)*max(1,norm(y));

                    if courbure > seuil
                        rhoBFGS=1/courbure;
                        B=(I-rhoBFGS*s*y')*B*(I-rhoBFGS*y*s') ...
                            +rhoBFGS*(s*s');
                        B=(B+B')/2;
                    else
                        B=I;
                    end
                end

                xval(:,i+1)=Xnew;
                result.crit(i+1)=Fnew;
                result.grad(i+1)=norm(Gnew);
                result.step(i)=alpha;
                result.iter=i;

                X=Xnew;
                F=Fnew;
                G=Gnew;

                if ~continuer
                    break
                end
            end

%%============ METHODE DE GAUSS-NEWTON ============
        case 'gauss-newton'
            for i=1:options.maxiter
                direction=resolutionStable(J'*J,-J'*r);

                if G'*direction >= 0 || any(~isfinite(direction))
                    direction=-G;
                end

                [pasAccepte,alpha,Xnew,Fnew,Gnew,~,rnew,Jnew]= ...
                    rechercheArmijo(X,F,G,direction,params,options);

                if ~pasAccepte
                    result.stop='RecherchePas';
                    result.iter=i-1;
                    break
                end

                [continuer,result]=testsortie( ...
                    X,Xnew,F,Fnew,Gnew,i,result,options);

                xval(:,i+1)=Xnew;
                result.crit(i+1)=Fnew;
                result.grad(i+1)=norm(Gnew);
                result.step(i)=alpha;
                result.iter=i;

                X=Xnew;
                F=Fnew;
                G=Gnew;
                r=rnew;
                J=Jnew;

                if ~continuer
                    break
                end
            end

%%============ METHODE DE LEVENBERG-MARQUARDT ============
        case 'levenberg-marquardt'
            mu=options.lmMuInitial;

            for i=1:options.maxiter
                pasAccepte=false;

                for tentative=1:options.maxBacktracking
                    muEssai=mu;
                    direction=resolutionStable( ...
                        J'*J+muEssai*eye(n),-J'*r);

                    if any(~isfinite(direction))
                        mu=min(2*mu,options.lmMuMax);
                        continue
                    end

                    Xnew=X+direction;
                    [Fnew,Gnew,~,rnew,Jnew]=rosenbrock(Xnew,params);

                    diminutionReelle=F-Fnew;
                    diminutionPredite=-(G'*direction ...
                        +direction'*(J'*J)*direction);

                    if diminutionPredite > 0 && isfinite(Fnew) ...
                            && all(isfinite(Gnew))
                        rhoLM=diminutionReelle/diminutionPredite;
                    else
                        rhoLM=-Inf;
                    end

                    if rhoLM > 0
                        pasAccepte=true;

                        if rhoLM >= 0.75
                            mu=max(mu/3,options.lmMuMin);
                        elseif rhoLM <= 0.25
                            mu=min(2*mu,options.lmMuMax);
                        end
                        break
                    end

                    mu=min(2*mu,options.lmMuMax);
                end

                if ~pasAccepte
                    result.stop='RecherchePas';
                    result.iter=i-1;
                    break
                end

                [continuer,result]=testsortie( ...
                    X,Xnew,F,Fnew,Gnew,i,result,options);

                xval(:,i+1)=Xnew;
                result.crit(i+1)=Fnew;
                result.grad(i+1)=norm(Gnew);
                result.step(i)=1;
                result.mu(i)=muEssai;
                result.iter=i;

                X=Xnew;
                F=Fnew;
                G=Gnew;
                r=rnew;
                J=Jnew;

                if ~continuer
                    break
                end
            end

        otherwise
            error('optimdescent:Methode','Methode inconnue: %s',options.method);
    end

    xh=X;
    result.temp=toc(tempsDebut);
end


function [pasAccepte,alpha,Xnew,Fnew,Gnew,Hnew,rnew,Jnew]= ...
        rechercheArmijo(X,F,G,direction,params,options)
%RECHERCHEARMIJO Rebroussement le long d'une direction fixee.

    pasAccepte=false;
    alpha=options.initialStep;
    deriveeDirectionnelle=G'*direction;

    Xnew=X;
    Fnew=F;
    Gnew=G;
    Hnew=[];
    rnew=[];
    Jnew=[];

    if ~isfinite(deriveeDirectionnelle) || deriveeDirectionnelle >= 0
        return
    end

    for rebroussement=1:options.maxBacktracking
        Xessai=X+alpha*direction;
        [Fessai,Gessai,Hessai,ressai,Jessai]=rosenbrock(Xessai,params);

        if isfinite(Fessai) && all(isfinite(Gessai)) ...
                && Fessai <= F+options.armijoC*alpha*deriveeDirectionnelle
            pasAccepte=true;
            Xnew=Xessai;
            Fnew=Fessai;
            Gnew=Gessai;
            Hnew=Hessai;
            rnew=ressai;
            Jnew=Jessai;
            return
        end

        alpha=options.backtrackFactor*alpha;
    end
end


function [continuer,result]=testsortie( ...
        X,Xnew,F,Fnew,Gnew,i,result,options)
%TESTSORTIE Teste les criteres uniquement sur un itere accepte.

    continuer=true;
    variationX=norm(Xnew-X)/max(1,norm(X));
    variationF=abs(Fnew-F)/max(1,abs(F));

    if norm(Gnew) <= options.tolG
        result.stop='TolG';
        continuer=false;
    elseif variationX <= options.tolX
        result.stop='TolX';
        continuer=false;
    elseif variationF <= options.tolF
        result.stop='TolF';
        continuer=false;
    end

    if ~continuer
        result.iter=i;
    end
end


function divergent=estDivergent(F,G,critereInitial,options)
    divergent=~isfinite(F) || any(~isfinite(G)) ...
        || abs(F) > options.divergenceFactor*critereInitial;
end


function direction=directionNewton(H,G)
%DIRECTIONNEWTON Regularise la Hessienne si elle n'est pas definie positive.

    if isempty(H) || any(~isfinite(H(:)))
        direction=-G;
        return
    end

    Hsym=(H+H')/2;
    [~,indefinie]=chol(Hsym);

    if indefinie == 0
        direction=-(Hsym\G);
    else
        decalage=max(1e-8,-min(eig(full(Hsym)))+1e-8);
        Hsym=Hsym+decalage*eye(length(G));
        direction=-(Hsym\G);
    end
end


function solution=resolutionStable(A,b)
%RESOLUTIONSTABLE Resout A*x=b apres une faible regularisation diagonale.

    A=(A+A')/2;
    decalage=1e-12*max(1,norm(A,'fro'));
    A=A+decalage*eye(size(A));
    solution=A\b;
end
