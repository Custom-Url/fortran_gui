


mat = cell(3,1);

% 1 - 2

position = [6 9 13 18] ;
zone = [1 1 1;2 1 1;3 2 1;4 2 1];
nzone = 4;
loi = 'voigt';
phases = [1 2];
fv = [0.5 0.5;0.5 0.5;0.5 0.5;0.5 0.5];
normale = [1 1 0;1 1 0;1 -1 0;1 -1 0];
direction = [0 0 1 ; 0 0 1;0 0 1 ; 0 0 1];

mat{1} = MateriauComp_FFTP(position,zone,nzone,loi,phases,fv,normale,direction) ;


% 2 - 2

position = [11 15 19] ;
zone = [1 1 2;2 1 2;3 1 2];
nzone = 3;
loi = 'voigt';
phases = [2 2];
fv = [3/8 5/8;5/8 3/8;7/8 1/8];
normale = [4 -1 0;4 -1 0;4 -1 0];
direction = [0 0 1;0 0 1;0 0 1];

mat{2} = MateriauComp_FFTP(position,zone,nzone,loi,phases,fv,normale,direction) ;


% 1 - 3

position = [4] ;
zone = [1 1 1];
nzone = 1;
loi = 'voigt';
phases = [1 3];
fv = [0.5 0.5];
normale = [1 1 0];
direction = [0 0 1];

mat{3} = MateriauComp_FFTP(position,zone,nzone,loi,phases,fv,normale,direction) ;


Rep = 'Composite2_mix';


Creation_fichiers_composites( mat,Rep ) ;
