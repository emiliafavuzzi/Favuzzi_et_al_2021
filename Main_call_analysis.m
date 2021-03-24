clear
disp('running KO2_brain2');
disp('slice1_side1');
[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('KO2_brain2/slice1_side1');
functiongenerateCountMatrixOnClusterMoreTmem('KO2_brain2/slice1_side1',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
fastintensity('KO2_brain2/slice1_side1',cellmaskTmem119);
clear
disp('slice1_side2');
[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('KO2_brain2/slice1_side2');
functiongenerateCountMatrixOnClusterMoreTmem('KO2_brain2/slice1_side2',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
 fastintensity('KO2_brain2/slice1_side2',cellmaskTmem119);

clear
disp('slice2_side1');

[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('KO2_brain2/slice2_side1');
fastintensity('KO2_brain2/slice2_side1',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
functiongenerateCountMatrixOnClusterMoreTmem('KO2_brain2/slice2_side1',cellmaskTmem119);

clear
disp('slice2_side2');

[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('KO2_brain2/slice2_side2');
functiongenerateCountMatrixOnClusterMoreTmem('KO2_brain2/slice2_side2',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
fastintensity('KO2_brain2/slice2_side2',cellmaskTmem119);

clear
disp('slice3_side1');

[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('KO2_brain2/slice3_side1');
functiongenerateCountMatrixOnClusterMoreTmem('KO2_brain2/slice3_side1',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
fastintensity('KO2_brain2/slice3_side1',cellmaskTmem119);
clear
disp('slice3_side2');

[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('control1_brain2/slice3_side2');
functiongenerateCountMatrixOnClusterMoreTmem('control1_brain2/slice3_side2',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
fastintensity('control1_brain2/slice3_side2',cellmaskTmem119);
clear
clear
disp('running KO1_brain3');
disp('slice1_side1');
[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('KO1_brain3/slice1_side1');
functiongenerateCountMatrixOnClusterMoreTmem('KO1_brain3/slice1_side1',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
fastintensity('KO1_brain3/slice1_side1',cellmaskTmem119);
clear
disp('slice1_side2');
[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('KO1_brain3/slice1_side2');
functiongenerateCountMatrixOnClusterMoreTmem('KO1_brain3/slice1_side2',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
 fastintensity('KO1_brain3/slice1_side2',cellmaskTmem119);
% 
clear
disp('slice2_side1');

[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('KO1_brain3/slice2_side1');
fastintensity('KO1_brain3/slice2_side1',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
functiongenerateCountMatrixOnClusterMoreTmem('KO1_brain3/slice2_side1',cellmaskTmem119);

clear
disp('slice1_side2');

[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('KO1_brain3/slice2_side2');
functiongenerateCountMatrixOnClusterMoreTmem('KO1_brain3/slice2_side2',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
fastintensity('KO1_brain3/slice2_side2',cellmaskTmem119);

clear
disp('slice3_side1');

[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('KO1_brain3/slice3_side1');
functiongenerateCountMatrixOnClusterMoreTmem('KO1_brain3/slice3_side1',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
fastintensity('KO1_brain3/slice3_side1',cellmaskTmem119);
clear
disp('slice3_side2');

[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('KO1_brain3/slice3_side2');
functiongenerateCountMatrixOnClusterMoreTmem('KO1_brain3/slice3_side2',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
fastintensity('KO1_brain3/slice3_side2',cellmaskTmem119);
clear
%

clear
disp('running control2_brain3');
disp('slice1_side1');
[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('control2_brain3/slice1_side1');
functiongenerateCountMatrixOnClusterMoreTmem('control2_brain3/slice1_side1',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
fastintensity('control2_brain3/slice1_side1',cellmaskTmem119);
clear
disp('slice1_side2');
[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('control2_brain3/slice1_side2');
functiongenerateCountMatrixOnClusterMoreTmem('control2_brain3/slice1_side2',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
 fastintensity('control2_brain3/slice1_side2',cellmaskTmem119);

clear
disp('slice2_side1');

[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('control2_brain3/slice2_side1');
fastintensity('control2_brain3/slice2_side1',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
functiongenerateCountMatrixOnClusterMoreTmem('control2_brain3/slice2_side1',cellmaskTmem119);

clear
disp('slice2_side2');

[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('control2_brain3/slice2_side2');
functiongenerateCountMatrixOnClusterMoreTmem('control2_brain3/slice2_side2',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
fastintensity('control2_brain3/slice2_side2',cellmaskTmem119);

clear
disp('slice3_side1');

[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('control2_brain3/slice3_side1');
functiongenerateCountMatrixOnClusterMoreTmem('control2_brain3/slice3_side1',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
fastintensity('control2_brain3/slice3_side1',cellmaskTmem119);
clear
disp('slice3_side2');

[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('control2_brain3/slice3_side2');
functiongenerateCountMatrixOnClusterMoreTmem('control2_brain3/slice3_side2',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
fastintensity('control2_brain3/slice3_side2',cellmaskTmem119);
clear
%clear

clear
disp('trying KO1_brain3 slice1_side2');
[cellmaskTmem119,cellmaskallTmem_Gabbr1_2]=functionmaskeMasksOnClusterMoreTmem('KO1_brain3/slice1_side2');
functiongenerateCountMatrixOnClusterMoreTmem('KO1_brain3/slice1_side2',cellmaskTmem119,cellmaskallTmem_Gabbr1_2);
 fastintensity('KO1_brain3/slice1_side2',cellmaskTmem119);
% 
clear
