This repository contains:
1) the source code for SBL-Thickness (in src/)
2) an example to calculate SBL thickness in sphere phantoms (in Sphere-Example/)

To use run the sphere example,
1) Clone this repository to your local computer
2) Download and compile the external dependencies. These include:
   - VTK (https://vtk.org/download/)
   - Armadillo (https://arma.sourceforge.net/download.html)
   - Voro++ (https://math.lbl.gov/voro++/download/)
3) Compile this repository. You may need to tell cmake where the VTK/arma/voro directories are.. I find this easiest to do by compiling with ccmake
4) Make sure that the paths listed in Sphere-Example/phantom-sphere.sh match the path where you cloned the repo.
5) Enter the Sphere-Example directory and run ./phantom-sphere
6) If you want to visualize the results, you can load the Sphere-Example/Surface-files/surfaceOuter_thickness.vtp into any program that can read a .vtp file. I like using 3D-Slicer (https://download.slicer.org/) or Paraview (https://www.paraview.org/download/)

Note:  we are including the binary executable with explicit permission, and it should not be redistributed for any purposes other than within the SBL thickness pipeline.
