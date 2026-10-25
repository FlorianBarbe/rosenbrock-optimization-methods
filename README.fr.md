# Méthodes d'optimisation sur la fonction de Rosenbrock

[English](README.md) | [Français](README.fr.md)

Ce projet propose une implémentation MATLAB, réalisée sans solveur
d'optimisation externe, de plusieurs méthodes classiques d'optimisation non
linéaire appliquées à la fonction de Rosenbrock généralisée. Il s'appuie sur le
premier exercice d'un sujet pédagogique de calcul scientifique, mais il est
présenté comme une étude numérique reproductible.

## Méthodes implémentées

- descente de gradient à pas fixe ;
- descente de gradient avec rebroussement d'Armijo ;
- gradient conjugué non linéaire de Fletcher-Reeves ;
- Newton à pas unitaire ;
- Newton sécurisé avec recherche d'Armijo ;
- quasi-Newton inverse-BFGS ;
- Gauss-Newton ;
- Levenberg-Marquardt avec amortissement adaptatif.

Les méthodes partagent les mêmes critères d'arrêt et le même format
d'historique. Le code distingue l'arrêt par norme du gradient, la stagnation,
la limite d'itérations, l'échec de recherche de pas et la divergence.

## Problème étudié

La fonction généralisée est

$$
f(x)=\sum_{i=1}^{n-1}\left[b(x_{i+1}-x_i^2)^2+(1-x_i)^2\right],
$$

avec $b>0$. Son unique minimum global est atteint en
$x^\star=(1,\ldots,1)^T$ et vaut zéro.

Les expériences utilisent $n=2$ et $b=2$. L'implémentation fournit également
le gradient, la Hessienne, le vecteur des résidus et sa Jacobienne, avec la
convention cohérente

$$
f(x)=r(x)^Tr(x), \qquad \nabla f(x)=2J(x)^Tr(x).
$$

## Prérequis

- MATLAB, testé avec MATLAB R2025b ;
- aucune dépendance à Optimization Toolbox ou Symbolic Math Toolbox ;

## Exécution

Depuis la racine du projet dans MATLAB :

```matlab
main
```

Ce script vérifie les dérivées par différences finies, produit les
visualisations, étudie les trois pas fixes et compare les sept configurations
depuis les deux points initiaux.

Pour régénérer séparément les figures et les tableaux CSV :

```matlab
export_figures
```

L'exécution normale de `main.m` n'exporte aucun fichier.

## Paramètres principaux

| Paramètre | Valeur |
|---|---:|
| Coefficient de Rosenbrock $b$ | 2 |
| Tolérance sur le gradient | $10^{-8}$ |
| Tolérance relative sur les itérés | $10^{-8}$ |
| Tolérance relative sur le critère | $10^{-8}$ |
| Nombre maximal d'itérations | 1000 |
| Constante d'Armijo | $10^{-4}$ |
| Facteur de rebroussement | 0,75 |
| Pas d'essai initial | 1 |
| Pas fixes | $1$, $10^{-2}$, $10^{-4}$ |

## Résultats essentiels

| Méthode | $x_0=(-2,10)^T$ | $x_0=(4,-5)^T$ |
|---|---:|---:|
| Gradient Armijo | 903 (`TolF`) | 149 (`TolF`) |
| Fletcher-Reeves | 57 (`TolF`) | 74 (`TolF`) |
| Newton à pas unitaire | 0 (`Divergence`) | 5 (`TolG`) |
| Newton Armijo | 13 (`TolG`) | 10 (`TolG`) |
| BFGS | 18 (`TolG`) | 22 (`TolF`) |
| Gauss-Newton | 3 (`TolG`) | 2 (`TolG`) |
| Levenberg-Marquardt | 5 (`TolF`) | 4 (`TolG`) |

Un arrêt `TolF` indique une faible variation du critère et n'est pas présenté
automatiquement comme une convergence du gradient. La divergence de Newton à
pas unitaire depuis le premier point est conservée dans les résultats.

![Comparaison depuis le premier point initial](figures/partie1_06_comparaison_rosenbrock_point_1.png)

![Comparaison depuis le second point initial](figures/partie1_07_comparaison_rosenbrock_point_2.png)

## Organisation

```text
.
├── main.m
├── optimdescent.m
├── rosenbrock.m
├── verify_derivatives.m
├── export_figures.m
├── figures/
├── rapport/Rapport.pdf
└── csopt_tp_2020.pdf
```

## Limites

- La comparaison se limite au cas bidimensionnel avec $b=2$ et à deux points
  initiaux.
- Les temps sont des mesures uniques et trop courts pour établir un classement
  fiable des performances.
- Les critères de stagnation peuvent arrêter une méthode avant que la tolérance
  sur le gradient ne soit satisfaite.
- Le code a une vocation pédagogique et ne remplace pas une bibliothèque
  d'optimisation industrielle.

Les démonstrations, tableaux complets et commentaires sont disponibles dans
[`rapport/Rapport.pdf`](rapport/Rapport.pdf). Le fichier `csopt_tp_2020.pdf` est
conservé comme source pédagogique ; seul son exercice 1 est traité ici.

