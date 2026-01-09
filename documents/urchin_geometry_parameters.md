# Urchin Geometry Parameters

This document outlines the configuration parameters used to define the geometry of the Urchin model. Parameters are listed in the order they are computed in the generation workflow, starting from user inputs to derived global settings, and finally per-spike geometric derivations.

## 1. Input Parameters

User-defined settings that drive the generation process.

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `rCore` | Input | `1` | **Core Radius (nm)**. Radius of the central spherical core. |
| `spikeLength` | Input | `1` | **Spike Length (nm)**. Nominal length from core surface. |
| `spikeCount` | Input | `100` | **Number of Spikes**. Total spikes to generate. |
| `spikeTip` | Input | `[]` | **Spike Tip Diameter (nm)**. Diameter of tip cap. Defaults to `spikeLength/10` if empty. |
| `spikeConicality` | Input | `0.5` | **Conicality [-1,1]**. Controls shape: <br> `1`: Widest cone (max base). <br> `0`: Cylindrical stud. <br> `-1`: Tapered base (radius 0). |
| `resolution` | Input | `100` | **Mesh Resolution**. Dimensionless factor controlling mesh density. |
| `useFillet` | Input | `true` | **Enable Fillet**. Enables toroidal blending between core and spikes. |
| `flucFactor` | Input | `0.5` | **Fluctuation Factor** `[0,1]`. Magnitude of length variation. |
| `flucMethod` | Input | `'uniform'` | **Fluctuation Method**. `uniform` (Sobol), `random`, or `gaussian`. |
| `distMethod` | Input | `'uniform'` | **Distribution Method**. `uniform` (Golden Spiral) or `random`. |
| `refinedOrientation`| Input | `true` | **Refine Orientation**. Enable Coulomb-like relaxation of positions. |
| `refinedOrientationThreshold`| Input | `0.1` | **Refinement Threshold**. Min angular separation target (degrees). |

## 2. Global Derived Parameters

Parameters derived immediately from inputs that affect the entire mesh.

| Parameter | Definition / Derivation | Description |
| :--- | :--- | :--- |
| `scale` | `2 * (rCore + spikeLength)` | **Characteristic Scale**. Normalization scale for mesh heuristics. |
| `minSpacing` | `max(1e-9, scale / resolution)` | **Minimum Spacing**. Target minimum vertex distance. |
| `coreSubdiv` | `ceil(0.5 * log2((Area / spacing^2))` | **Core Subdivision Level**. Icosphere recursion depth based on resolution. |

## 3. Spike Orientations and Lengths

First stage of per-spike computation.

| Parameter | Definition / Derivation | Description |
| :--- | :--- | :--- |
| `spikeOrientations` | `distMethod` + Relaxation | **Spike Directions**. Unit vectors $\vec{n}_i$ for spike axes. |
| `spikeLengths` | `spikeLength` + `flucFactor` | **Actual Spike Lengths**. $L_i$. Final length per spike after fluctuations. |
| `zTipApex` | `rCore + spikeLengths` | **Tip Apex Distance**. Distance from origin to the very top of the spike. |

## 4. Spike Base Constraints

Limits calculated to prevent spike overlap and ensure geometric validity.

| Parameter | Definition / Derivation | Description |
| :--- | :--- | :--- |
| `thetaMins` | Nearest Neighbor Dot Product | **Nearest Neighbor Angle**. Smallest angle to any other spike. |
| `alphaBaseMaxSpike` | `0.5 * thetaMins` | **Max Base Angle (Neighbors)**. Limit to prevent overlap with neighbors. |
| `alphaBaseMaxTangent` | `acos((rCore - rTip) / (zApex - rTip))` | **Max Base Angle (Tangency)**. Limit for core tangency (horizon). |
| `alphaBaseMax` | `min(alphaBaseMaxSpike, ...)` | **Effective Max Base Angle**. The stricter of neighbor or tangent limits. |
| `rBaseMax` | `rCore * sin(alphaBaseMax)` | **Maximum Base Radius**. Radius corresponding to max angle. |

## 5. Spike Tip Geometry

Defining the spherical cap at the tip.

| Parameter | Definition / Derivation | Description |
| :--- | :--- | :--- |
| `rTipMax` | `zTipApex / (1 + 1/sin(alphaBaseMax))` | **Maximum Tip Radius**. Largest tip radius possible within the max cone angle. |
| `rTip` | `min(inputSpikeTip/2, rTipMax)` | **Actual Tip Radius**. Input radius clamped by geometric limits. |
| `zTipCenter` | `zTipApex - rTip` | **Tip Sphere Center**. Distance to center of tip sphere. |

## 6. Spike Base Geometry

Defining the connection to the core, controlled by `conicality`.

| Parameter | Definition / Derivation | Description |
| :--- | :--- | :--- |
| `rBase` | *Conicality Interpolation* | **Actual Base Radius**. <br> If `C >= 0`: `rTip + C * (rBaseMax - rTip)` <br> If `C < 0`: `(1 + C) * rTip` |
| `alphaBase` | `asin(rBase / rCore)` | **Base Angle**. Polar angle of the spike base on core. |
| `zBase` | `sqrt(rCore^2 - rBase^2)` | **Base Height**. Distance from origin to the base plane. |

## 7. Seam Geometry (Cone-Tip Tangency)

Solving for the tangential connection between the conical body and spherical tip.

| Parameter | Definition / Derivation | Description |
| :--- | :--- | :--- |
| `zTipSeam` | Derived via System Solver | **Seam Height**. Z-coord where cone meets tip sphere tangentially. |
| `rTipSeam` | `sqrt(rTip^2 - (zTipSeam - zTipCenter)^2)` | **Seam Radius**. Radius of the ring at the seam. |
| `alphaCone` | `atan((rBase - rTipSeam)/(zTipSeam - zBase))` | **Cone Half-Angle**. Slope of the conical section. |

### Seam Derivation System
The seam parameters (`zTipSeam`, `rTipSeam`, `alphaCone`) are the solution to the system of equations ensuring tangency:
1. `tan(alphaCone) = (rBase - rTipSeam) / (zTipSeam - zBase)` (Cone Slope)
2. `tan(alphaCone) = (zTipSeam - zTipCenter) / rTipSeam` (Sphere Tangent)
3. `rTip^2 = (zTipSeam - zTipCenter)^2 + rTipSeam^2` (Sphere Surface)
