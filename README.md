# Urchin Creator

**Urchin Creator** is a MATLAB-based tool for generating 3D models of spherical urchin-like nanoparticles with customizable conical spikes. These models are suitable for use in scientific simulations, such as electromagnetic or fluid dynamics studies.

The core of this project is the `urchin.m` function, which allows for detailed control over the urchin's geometry, including core size, spike length, number of spikes, tip thickness, conicality, and surface smoothing.

A Python port of this tool is planned for the future.

## Project Structure

-   `matlab/src/`: Contains the main MATLAB source code (`urchin.m`).
-   `matlab/examples/`: Includes example scripts demonstrating how to use `urchin.m`.
-   `matlab/tests/`: Contains unit tests for the MATLAB code.
-   `python/`: Placeholder for the future Python version.
-   `docs/`: Project documentation.
-   `data/`: Directory for example output files.
-   `.github/workflows/`: GitHub Actions for Continuous Integration.

## MATLAB Version

### Dependencies
-   MATLAB
-   Image Processing Toolbox
-   Statistics and Machine Learning Toolbox (for `sobolset`, `net` if using 'uniform' `flucMethod`)
-   (Optional) Parallel Computing Toolbox (for GPU acceleration)
-   (For visualization in examples/main function) `viewer3d`, `volshow`, `surfaceMeshShow` (which might be part of a specific toolbox or require separate installation/availability).

### How to Use
1.  Ensure MATLAB and the required toolboxes are installed.
2.  Add the `matlab/src` directory to your MATLAB path.
    ```matlab
    addpath('path/to/UrchinCreator/matlab/src');
    ```
3.  Run the example script `matlab/examples/run_urchin_creation.m` or call the `urchin` function directly:
    ```matlab
    % Example call
    [mesh, mask, threshold, eqRadius] = urchin('cr', 30, 'sl', 15, 'ns', 100, 'st', 2);
    
    % To visualize if no output arguments are requested
    urchin('cr', 30, 'sl', 15, 'ns', 100, 'st', 2);
    ```

Refer to the documentation in `docs/` and the comments within `matlab/src/urchin.m` for detailed information on parameters.

## Citation
If you use this software in your research, please cite it using the information in `CITATION.cff`.

## License
This project is licensed under the MIT License - see the `LICENSE` file for details.
