C***********************************
C DEFINITION DES SYST DE GLISST
C de referencee
c avant rotation dans les grains
c valables pour tous les grains 
c     definir le nb de syst Nsys et les coord
c     de xmsys(i) et xnsys(i) pour le syst i
c     en data dans la procedure
c les syst sont donnes dans le repere cristallin
c
c ordre : les n dans XNSYS puis les m dans XMSYS
c    d'abord les coord selon e1, puis celles selon e2 et enfin celles selon e3
c
c*****************************************************************************
C
      DATA XNSYS 
c 1ere coordonnee
     1            /0.70710678,-0.70710678,0.,
     2            0.,0.70710678,-0.70710678,
     3           -0.70710678,0.70710678,0.,
     4           0.70710678,-0.70710678,0.,
c 2ieme coordonnee
     7           0.70710678,0.,-0.70710678,
     8           -0.70710678,0.,0.70710678,
     9           -0.70710678,0.,0.70710678,
     1           -0.70710678,0.,0.70710678,
c 3ieme coordonnee
     4          0.,0.70710678,-0.70710678,
     5          0.70710678,-0.70710678,0.,
     6          0.,0.70710678,-0.70710678,
     7          0.,-0.70710678,0.70710678/
c
      DATA XMSYS
c 1ere coordonnee
     1           /0.57735027,0.57735027,0.57735027,
     2          0.57735027,0.57735027,0.57735027,
     3          -0.57735027,-0.57735027,-0.57735027,
     4          0.57735027,0.57735027,0.57735027,
c 2ieme coordonnee
     7          -0.57735027,-0.57735027,-0.57735027,
     8          0.57735027,0.57735027,0.57735027,
     9          0.57735027,0.57735027,0.57735027,
     1          0.57735027,0.57735027,0.57735027,
c 3ieme coordonnee
     4          0.57735027,0.57735027,0.57735027,
     5          0.57735027,0.57735027,0.57735027,
     6          0.57735027,0.57735027,0.57735027,
     7          -0.57735027,-0.57735027,-0.57735027/
c
