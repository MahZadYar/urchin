## urchin - MATLAB Example

% This script demonstrates how to use the urchin.m function to generate
% a 3D urchin model with specified parameters and export it.

% Add the source directory to MATLAB path (if not already)
% addpath('../src'); % Assuming this script is in urchin/matlab/examples/

% [text] INPUT PARAMETERS
cr    = 32.5;  % Core Radius (nm)
hr    = 0;     % Hollow Radius (nm) - set to 0 for a solid core
sl    = 15;    % Spikes Length (nm)
ns    = 200;   % Number of Spikes
st    = 5;     % Spikes tip thickness (nm)
sc    = 0.5;   % Conicality of spikes (-1 to 1)
sf    = 0.5;   % Spikes length fluctuation (0 to 1)
flucMethod = 'uniform'; % 'uniform' or 'random'
distMethod = 'uniform'; % 'uniform' or 'random'
smth  = 1;     % Smoothing factor
res   = 0.5;   % Voxel dimension (nm)
antialiasing = true; % Enable anti-aliasing

% [text] FILE NAMING & OUTPUT PATH
% Define a base folder for outputs, e.g., within the project or a specific data folder
% For this example, we'll suggest creating an 'output' folder within 'matlab/examples'
outputBaseFolder = fullfile(pwd, 'output'); % Creates 'output' in the current example script directory
if ~exist(outputBaseFolder, 'dir')
    mkdir(outputBaseFolder);
end

filenamePattern = sprintf('Urchin_cr%g_sl%g_ns%d_st%g_sc%g_sf%g_res%g', cr, sl, ns, st, sc, sf, res);
filename = matlab.lang.makeValidName(filenamePattern); % Ensure filename is valid

% [text] OPTIONS FOR VISUALIZATION AND EXPORT
vis_mask   = false; % Visualize the urchin mask (can be memory intensive for high res)
exp_mask   = true;  % Export Mask as .mat file
vis_mesh   = true;  % Visualize the urchin mesh
exp_stl    = true;  % Export Mesh as .stl file

fprintf('Starting Urchin Generation with the following parameters:\n');
fprintf('  Core Radius (cr): %.2f nm\n', cr);
fprintf('  Hollow Radius (hr): %.2f nm\n', hr);
fprintf('  Spike Length (sl): %.2f nm\n', sl);
fprintf('  Number of Spikes (ns): %d\n', ns);
fprintf('  Spike Tip Thickness (st): %.2f nm\n', st);
fprintf('  Conicality (sc): %.2f\n', sc);
fprintf('  Spike Fluctuation (sf): %.2f\n', sf);
fprintf('  Fluctuation Method: %s\n', flucMethod);
fprintf('  Distribution Method: %s\n', distMethod);
fprintf('  Smoothing (smth): %.2f\n', smth);
fprintf('  Resolution (res): %.2f nm\n', res);
fprintf('  Antialiasing: %s\n', mat2str(antialiasing));
fprintf('Output filename prefix: %s\n', filename);
fprintf('Output folder: %s\n', outputBaseFolder);

% [text] BUILD URCHIN
try
    if nargout(@urchin) == 0 && (vis_mask || vis_mesh) % If urchin handles its own visualization
        fprintf('\nCalling urchin function for direct visualization...\n');
        urchin('cr', cr, 'hr', hr, 'sl' , sl ,'ns', ns ,'st', st ,'sc', sc ,'sf', sf ,'res', res ,'smth', smth ,'flucMethod', flucMethod, 'distMethod', distMethod, 'antialiasing', antialiasing);
        if exp_mask || exp_stl
            fprintf('Warning: Direct visualization call. Re-run with output arguments to export.\n')
        end
    else
        fprintf('\nCalling urchin function to get outputs...\n');
        [mesh, mask, threshold, eqRadius] = urchin('cr', cr, 'hr', hr, 'sl' , sl ,'ns', ns ,'st', st ,'sc', sc ,'sf', sf ,'res', res ,'smth', smth ,'flucMethod', flucMethod, 'distMethod', distMethod, 'antialiasing', antialiasing);
        fprintf('Urchin generation complete. Equivalent Radius: %.2f nm, Threshold: %.4f\n', eqRadius, threshold);

        % [text] VISUALIZATION (if not handled by urchin directly)
        if vis_mask
            fprintf('Visualizing mask...\n');
            figure;
            volshow(mask, 'Colormap', parula, 'RenderingStyle', 'CinematicRendering');
            title(sprintf('Urchin Mask - %s', strrep(filename, '_', ' ')));
        end

        if vis_mesh
            fprintf('Visualizing mesh...\n');
            figure;
            if exist('surfaceMeshShow', 'file')
                surfaceMeshShow(mesh, 'Wireframe', false, 'FaceColor', [0.8 0.7 0.6]);
            elseif isfield(mesh, 'Vertices') && isfield(mesh, 'Faces')
                patch('Vertices', mesh.Vertices, 'Faces', mesh.Faces, 'FaceColor', [0.8 0.7 0.6], 'EdgeColor', 'none', 'FaceLighting', 'gouraud');
                axis equal; view(3); camlight; lighting gouraud;
            else
                warning('Could not visualize mesh. surfaceMeshShow not found and mesh structure unrecognized.');
            end
            title(sprintf('Urchin Mesh - %s', strrep(filename, '_', ' ')));
        end
        
        % [text] EXPORT MASK
        if exp_mask
            mask_filename = fullfile(outputBaseFolder, [filename, '_mask.mat']);
            fprintf('Exporting mask to: %s\n', mask_filename);
            save(mask_filename, 'mask', 'cr', 'sl', 'ns', 'st', 'sc', 'sf', 'res', 'smth', 'flucMethod', 'distMethod', 'antialiasing', 'hr', 'eqRadius', 'threshold');
        end

        % [text] EXPORT MESH
        if exp_stl
            stl_filename = fullfile(outputBaseFolder, [filename, '_mesh.stl']);
            fprintf('Exporting STL mesh to: %s\n', stl_filename);
            if exist('writeSurfaceMesh', 'file')
                writeSurfaceMesh(mesh, stl_filename, 'Encoding', 'binary');
            elseif isfield(mesh, 'Vertices') && isfield(mesh, 'Faces') && exist('stlwrite', 'file') % Check for stlwrite from File Exchange
                stlwrite(stl_filename, mesh.Faces, mesh.Vertices);
            elseif isfield(mesh, 'Vertices') && isfield(mesh, 'Faces') % Basic stlwrite if available
                 try 
                    TR = triangulation(double(mesh.Faces), double(mesh.Vertices));
                    stlwrite(TR, stl_filename, 'binary');
                 catch ME
                    warning('Failed to export STL. writeSurfaceMesh not found, stlwrite (File Exchange) not found or failed. Error: %s', ME.message);
                 end
            else
                warning('Could not export STL. writeSurfaceMesh not found and mesh structure unrecognized or stlwrite not available.');
            end
        end
    end
    fprintf('\nExample script finished.\n');
catch ME
    fprintf(2, 'An error occurred: %s\n', ME.message);
    fprintf(2, 'Error in %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
end
