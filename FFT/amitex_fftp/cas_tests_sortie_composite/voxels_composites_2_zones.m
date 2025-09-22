


mat = cell(1,1);

% 1 - 1

position = [6 9 11 15 19 13 18] ;
zone = [1 1 3;2 1 3;3 2 3;4 2 3;5 2 3;6 3 4;7 3 4];
nzone = 7;
loi = 'voigt';
phases = [1 1];
fv = [0.5 0.5;0.5 0.5;5/8 3/8 ; 3/8 5/8 ; 1/8 7/8;0.5 0.5;0.5 0.5];
normale = [1 1 0;1 1 0;4 -1 0;4 -1 0;4 -1 0;1 -1 0;1 -1 0];
direction = [0 0 1 ; 0 0 1;0 0 1;0 0 1;0 0 1;0 0 1 ; 0 0 1];

mat{1} = MateriauComp_FFTP(position,zone,nzone,loi,phases,fv,normale,direction) ;



Rep = 'Composite2zones';


Creation_fichiers_composites( mat,Rep ) ;
