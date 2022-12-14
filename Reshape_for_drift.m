%New rearanging of 3D array in order to account for drift with each loop. 

Array1
Array2

% Upload both arrays
% Once uploaded make Array2 equal to Array 1
Array2 = Array2(newsizex,newsizey,newsizez);
% Drift correct 
% these are the pixel values 

True_drift_x
True_drift_y
True_drift_z

% Modify the Array1 by drift amount
Array1 = Array1(1:(xsize-True_drift_x),1:(ysize - True_drift_y),1:(zsize-True_drift_z));

Array1 = Array1(1:(xsize - 10), 1:(ysize - 10), 1:(zsize - 10));
% modify Array2 on opposite side by Drift amount. 
Array2 = Array2((1 + True_drift_x) : xsize, (1 + True_drift_y) : ysize, (1 + True_drift_z) : zsize);

% Perform 3D PIV.

[newsizex,newsizey,newsizez] = size(Array1); %Newsizes