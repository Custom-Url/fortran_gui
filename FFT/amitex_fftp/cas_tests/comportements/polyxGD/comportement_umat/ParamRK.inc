C=======================================================================
C=  Parametres utiles pour les sousprogrammes d'integration numerique  =
C=  d'un modele de comportement par une methode de RUNGE-KUTTA d'ordre =
C=  p avec determination du pas de temps automatique avec estimation   =
C=  de l'erreur par une methode de RUNGE-KUTTA d'ordre q (=p-1)        =
C=======================================================================

C= Nombre maximal de grandeurs pouvant etre integrees dans MPNRKpq
      PARAMETER ( nVloc = 81 )
C= Nombre maximal de parametres externes utilises dans MPNRKpq
      PARAMETER ( nPext = 30 )

C= Criteres de precision et de tolerance pour l'integration numerique
C=  - toly : tolerance sur les grandeurs (valeur absolue)
C=  - eps  : erreur (locale) d'integration recherchee
      PARAMETER ( toly = 1.D-4 , eps = 1.D-5 )

C= Nombre maximal de (sous)pas calcules pour l'integration sur un pas
      PARAMETER ( IterMax = 50000 )

C= Coefficient de securite utile au calcul de la taille d'un souspas
      PARAMETER ( hoptsec = 0.90 )

C= Valeurs min. et max. pour la taille du souspas suivant
      PARAMETER ( hoptmin = 0.02 , hoptmax = 10. )

C= Valeur pour le decoupage du souspas suite a erreur 999 lors du calcul
C= des derivees en un point de la methode RKpq
      PARAMETER ( hopt999 = 0.01 )

C= Exposant de l'erreur dans le cas de la methode RK21
      PARAMETER ( hc21pui = 0.50 )
C= Exposant de l'erreur dans le cas de la methode RK43
      PARAMETER ( hc43pui = 0.25 )

C= Pas du gradient de F au dela duquel on redecoupe en deux le pas
      PARAMETER (pasF = 1.D-4)