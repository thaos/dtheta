library(dtheta)
library(R.matlab)
library(future)

future::plan(strategy = "multicore")

## Systeme de Pomeau-Manneville avec les parametres par defaut du papier de Davide
"pom.man" = function(xt,a=0.91,epsil=0)
{
    xt=max(xt,0)
    X = (xt + xt^(1+a)) %% 1 + epsil*rnorm(1,0,1)
    return(X)
}

## Calcul d'une trajectoire de P-M avec une condition initale X0
"pom.man.traj" = function(X0,N=100,a=0.91,epsil=0)
{
    Xt=X0
    X=X0
    for(i in 1:(N-1)){
        Xt= pom.man(Xt,a,epsil)
        X=c(X,Xt)
    }
    invisible(X)
}

## Creation d'une trajectoire longue de reference, pour les calcul de d/theta
X0=0.1
N=10^5
X.ref=pom.man.traj(X0,N)
X.0=seq(0.001,1,by=0.001)

if(FALSE){
  ## Lecture des donnees de d/theta de Davide
  ## La trajectoire de reference est calculee avec la meme condition initiale X0=0.1
  ## Les trajectoires dtheta.df$xref et X.ref divergent au bout de 150 pas de temps,
  ## a cause d'erreurs d'arrondis
  filin = "D_theta_PM_Manneville_ref.mat"
  dtheta.df = R.matlab::readMat(filin)
  ## dtheta.df:
  ## > names(dtheta.df)
  ##[1] "D"      "a"      "quanti" "theta"  "xref"

  ## Calculs de d/theta sur les 1000 premiers points de X.ref en utilisant
  ## les points suivants
  ## Calculs avec la trajectoire de reference X.ref faite en R
  dX.Ferro.Xref = dtheta_allpoints(X.ref[1:1000],X.ref[1001:100000],method="Ferro")
  dX.Sueve.Xref = dtheta_allpoints(X.ref[1:1000],X.ref[1001:100000],method="Sueveges")

  ## Calculs avec la trajectoire de reference dtheta.df$xref de Davide, faite en matlab
  dX.Ferro.df = dtheta_allpoints(X.ref[1:1000],dtheta.df$xref[1,1001:100000],method="Ferro")
  dX.Sueve.df = dtheta_allpoints(X.ref[1:1000],dtheta.df$xref[1,1001:100000],method="Sueveges")

  ## Calculs sur X.0 avec la trajectoire de reference X.ref faite en R
  dX.0.Ferro.Xref = dtheta_allpoints(X.0,X.ref,method="Ferro")
  dX.0.Sueve.Xref = dtheta_allpoints(X.0,X.ref,method="Sueveges")

  ## Calculs sur X.0 avec la trajectoire de reference dtheta.df$xref de Davide, faite en matlab
  dX.0.Ferro.df = dtheta_allpoints(X.0,dtheta.df$xref[1,],method="Ferro")
  dX.0.Sueve.df = dtheta_allpoints(X.0,dtheta.df$xref[1,],method="Sueveges")
}
