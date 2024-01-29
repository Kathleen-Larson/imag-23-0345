#!/bin/bash
#set -x

PHANTOM="sphere"
WDIR=/home/larsonke/Code/SBL-Thickness/Sphere-Example ## This should match the path to the cloned repository


# Executables (compiled from src)
EXE_DIR=$WDIR/../bin
decimateMesh=${EXE_DIR}/decimateMesh
getMeshEdgeLengths=${EXE_DIR}/getMeshEdgeLengths
integrateFieldLines=${EXE_DIR}/integrateFieldLines
laplacianSolver=${EXE_DIR}/laplacianSolver
makeAnnulus=${EXE_DIR}/makeAnnulus
spmesh=${EXE_DIR}/spmesh
spmesh2voroInput=${EXE_DIR}/spmesh2voroInput
surface2spmeshInput=${EXE_DIR}/surface2spmeshInput
voro= ##<path to locally compiled voro++ executable>

# Set up directories
ANNULUS_DIR=$WDIR/Annulus-Mesh-Files
SPMESH_DIR=$WDIR/SPMESH-Files
SURFACE_DIR=$WDIR/Surface-Files
mkdir -p $WDIR $ANNULUS_DIR $SURFACE_DIR $SPMESH_DIR


########################
surfaceOuter_base=${PHANTOM}_outer
surfaceInner_base=${PHANTOM}_inner

surfaceOuter=${SURFACE_DIR}/${surfaceOuter_base}.vtp
surfaceInner=${SURFACE_DIR}/${surfaceInner_base}.vtp

if [ ! -f $surfaceOuter ] || [ ! -f $surfaceInner ] ; then
    "Cannot find input surfaces! exiting..."
    exit
fi


# Create annulus
echo "Building input annulus..."
annulus_base=annulus
annulus=$ANNULUS_DIR/${annulus_base}.vtp
$makeAnnulus $surfaceOuter $surfaceInner $annulus


# Convert polydata to spmesh input
spmeshInput=${SPMESH_DIR}/spmeshInput.nml
$surface2spmeshInput $annulus $spmeshInput "Phantom example (${PHANTOM})"


# Generate user input file for spmesh
spmeshOutput_base=${SPMESH_DIR}/spmeshOutput
spmeshUserInputFile=${SPMESH_DIR}/spmeshUserInput.txt
input_del=$($getMeshEdgeLengths $annulus "Median")

echo -n > $spmeshUserInputFile
echo 1 >> $spmeshUserInputFile #Input internal model size (1 - 10000)
echo 0 >> $spmeshUserInputFile #Input the max. number of auto-refinement level
echo $spmeshInput >> $spmeshUserInputFile #Input NML+ File Name
echo 8 >> $spmeshUserInputFile #Input Element Type (cube building block=4, deltahedral block=8)
echo 1 >> $spmeshUserInputFile #Pattern (0= Cube / 1= Deltahedraron)
echo $input_del >> $spmeshUserInputFile #Input Delta Size
echo 0 >> $spmeshUserInputFile #Need Refinement? (No: 0 / Yes: 1)
echo 0 >> $spmeshUserInputFile #Export Numerical Surfaces? (Yes=1/No=0)
echo 1 >> $spmeshUserInputFile #Need output ? (1=yes / 0=No)
echo 1 >> $spmeshUserInputFile #Output file format (WPI: 0/ Thayer: 1)
echo 0 >> $spmeshUserInputFile #Elem. numbering convension (normals in=0, normals out=1)
echo 0 >> $spmeshUserInputFile #File Type (Ascii: 0 | Binary: 1)
echo $spmeshOutput_base >> $spmeshUserInputFile #Output file name w/o extension
echo -1 >> $spmeshUserInputFile #How many boundary surfaces need be extracted [1-1 (-1 = quit)]


# Build FDM w/ spmesh
echo "Generating FDM..."
$spmesh < $spmeshUserInputFile >> /dev/null
mv sp* $SPMESH_DIR/


# Calculate Voronoi diagram of FDM
spmeshOutputUGrid=${spmeshOutput_base}.vtk
voronoiDiagram_base=${ANNULUS_DIR}/voronoi_diagram

bounds_list=$($spmesh2voroInput $spmeshOutput_base $annulus $surfaceOuter $surfaceInner $voronoiDiagram_base $spmeshOutputUGrid)
$voro -o -v -c %i" "%s" "%n" "%f $bounds_list $voronoiDiagram_base >> /dev/null


# Solve for Laplacian gradient over the ugrid
echo "Calculating Laplacian gradient..."
FDMugrid=${ANNULUS_DIR}/FDM.vtk

$laplacianSolver $spmeshOutputUGrid $annulus $surfaceOuter $surfaceInner ${voronoiDiagram_base}.vol $FDMugrid


# Integrate thickness
echo "Integrating thickness..."
surfaceOuter_thickness=$SURFACE_DIR/${surfaceOuter_base}_thickness.vtp
step_size=0.5
direction=1 # integration from outer-->inner (-1 for inner-->outer)

$integrateFieldLines $FDMugrid $surfaceOuter $surfaceInner $surfaceOuter_thickness $step_size $direction
