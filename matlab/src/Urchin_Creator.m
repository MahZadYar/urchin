% [text] INPUT PARAMETERS
cr    = 65/2;  % Core Radius (nm) %[control:editfield:24b3]{"position":[15,17]}
hr    = 0;  % Core Height (nm) %[control:editfield:24b3]{"position":[15,17]} 29.3
sl    = 15;   % Spikes Length (nm) %[control:editfield:2b3b]{"position":[15,17]}
ns    = 200;  % Number of Spikes %[control:spinner:7090]{"position":[15,17]}
st    = 5;   % Spikes tip thickness (0<=st<=rs) %[control:editfield:517e]{"position":[15,16]}
sc    = 0.5;   % Conicality of spikes -1(reversed conical)<=sc<=1(conical) %[control:slider:3f12]{"position":[15,16]}
sf    = 1; % Spikes fluctuation 0(uniform)<=sf<=1(maximally fluctuated) %[control:slider:3f12]{"position":[15,16]}
flucMethod = 'uniform'; % 'uniform', 'random'
distMethod = 'uniform'; % 'uniform', 'random'
smth  = 1; %[control:editfield:9196]{"position":[15,18]}
res = 0.5; % Voxel dimension (nm) %[control:editfield:34c2]{"position":[15,18]}
antialiasing = true; % Anti-aliasing (nm) %[control:editfield:34c2]{"position":[15,18]}

vis_mask   = false; % Visualize the urchin mask %[control:checkbox:5ea4]{"position":[20,25]}
exp_mask   = false; % Export Mask %[control:checkbox:6629]{"position":[20,25]}
vis_mesh   = true;  % Visualize the urchin mesh %[control:checkbox:80d6]{"position":[20,24]}
exp_stl    = true; % Export Mesh %[control:checkbox:9392]{"position":[19,24]}


% [text] FILE NAMING
filename = sprintf('Nano-Urchin_dc-%g_sl-%g_ns-%g-%s_st-%g_sc-%g_sf-%g-%s_hr-%g_res-%g_%s-%s', cr*2, sl, ns, distMethod, st, sc, sf, flucMethod, hr, res, smth * 'smth', antialiasing * 'aa');
foldername = 'D:\OneDrive - Kaunas University of Technology\~Science Projects\NanoTRAACES\Nano-Urchins Studies\Data\Simulations\Urchin Models';
if exp_mask | exp_stl
    % [text] BUILD URCHIN
    [mesh, mask, eqRadius] = urchin('cr', cr, 'hr', hr, 'sl' , sl ,'ns', ns ,'st', st ,'sc', sc ,'sf', sf ,'res', res ,'smth', smth ,'flucMethod', flucMethod, 'distMethod', distMethod, 'antialiasing', antialiasing);
    % [text] VISUALIZATION
    if vis_mask
        % colormap = [1 * linspace(0,1,256)', 0.75 * linspace(0,1,256)', 0 * linspace(0,1,256)'];
        volshow(mask_smoothed, colormap=parula, RenderingStyle="CinematicRendering");
    end

    if vis_mesh %[output:group:46c8f5f9]
        surfaceMeshShow(mesh, Wireframe=false); %[output:20a5132e]
    end %[output:group:46c8f5f9]
        
    % [text] EXPORT MASK
    if exp_mask
        save(fullfile(foldername, [filename, '.mat']), 'mask');
        % h5create(fullfile(foldername, [filename, '.h5']), '/mask', size(mask));
        % h5write(fullfile(foldername, [filename, '.h5']), '/mask', unit8(mask));
        % export as 3d density map?
        
    end

    if exp_stl
        writeSurfaceMesh(mesh, fullfile(foldername, [filename, '.stl']), "Encoding", "binary");
    end
else
    % [text] BUILD URCHIN
    urchin('cr', cr, 'hr', hr, 'sl' , sl ,'ns', ns ,'st', st ,'sc', sc ,'sf', sf ,'res', res ,'smth', smth ,'flucMethod', flucMethod, 'distMethod', distMethod, 'antialiasing', antialiasing);
end

