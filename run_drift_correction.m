

function [True_drift_x, True_drift_y, True_drift_z, Array_source, Array_search] = run_drift_correction (xsize, ysize, zsize, Array1, Array2)

%Center Of 3D Grid.
        centrex = ceil (xsize/2);
        centrey  = ceil (ysize/2);
        centrez = ceil (zsize/2);
        
        % Source size from center out for drift correction.
        x_source_start = ceil (centrex - (xsize*0.1));
        y_source_start = ceil (centrey - (ysize*0.1));
        z_source_start = ceil (centrez - (zsize*0.1));
        x_source_finish = ceil (centrex + (xsize*0.1));
        y_source_finish = ceil (centrey + (ysize*0.1));
        z_source_finish = ceil (centrez + (zsize*0.1));
        
        % Source size from center out for drift correction. Must be bigger than Search size.
        x_search_start = ceil (centrex - (xsize*0.2));
        y_search_start = ceil (centrey - (ysize*0.2));
        z_search_start = ceil (centrez - (zsize*0.2));
        x_search_finish = ceil (centrex + (xsize*0.2));
        y_search_finish = ceil (centrey + (ysize*0.2));
        z_search_finish = ceil (centrez + (zsize*0.2));
        disp("Run Drift Corr")
        % Building 3D Array from Source and search Size.
        Array_source = Array1 (x_source_start:x_source_finish, y_source_start:y_source_finish, z_source_start:z_source_finish);
        Array_search = Array2 (x_search_start:x_search_finish,y_search_start:y_search_finish,z_search_start:z_search_finish);
        disp("Drift Arraysource/search")
        % Convolution analysis between Source and Search Arrays for Drift
        % Corrction.
        drift_t_3d=convn(Array_source,Array_search(end:-1:1,end:-1:1,end:-1:1));
        disp("Drift Convo finished")
        % Offset
        xsize_drift_offset = size(Array_source,1) + size(Array_search,1);
        ysize_drift_offset = size(Array_source,2) + size(Array_search,2);
        zsize_drift_offset = size(Array_source,3) + size(Array_search,3);
        disp("Drift offset")
        % Finding Max of Convolution for Final Drift Value in XYZ.
        [drift_val_z,izz]= max(drift_t_3d,[],3);
        [drift_val_y,iyy]= max(drift_val_z,[],2);
        [drift_val_x,ixx]= max(drift_val_y,[],1);
        True_drift_x = ixx - (xsize_drift_offset/2);
        
        [drift_val_z,izz]= max(drift_t_3d,[],3);
        [drift_val_x,ixx]= max(drift_val_z,[],1);
        [drift_val_y,iyy]= max(drift_val_x,[],2);
        True_drift_y = iyy - (ysize_drift_offset/2);
        
        [drift_val_x,ixx]= max(drift_t_3d,[],1);
        [drift_val_y,iyy]= max(drift_val_x,[],2);
        [drift_val_z,izz]= max(drift_val_y,[],3);
        True_drift_z = izz - (zsize_drift_offset/2);
        disp("Run Drift Corr over")
end