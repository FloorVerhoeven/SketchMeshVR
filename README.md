# Project Title

Master thesis on a sketch-based 3D modeling system for VR.

## Getting Started
This git repository contains both a VR (VR-SketchMesh) and non-VR (SketchMesh) version of the sketch-based 3D modeling software. SketchMesh depends on the included version of libigl and VR-SketchMesh depends on the specially adapted version of libigl called libiglVR plus the OculusSDK, OVRAvatarSDK and OVRPlatformSDK. 

### Prerequisites

The VR version is tested on Windows 10 with Microsoft Visual Studio 15 2017 (with Visual C++), the non-VR version is additionally also tested on OS X El Capitan with clang 8.0.0.

The Oculus Setup software should be installed on the PC.

### Installing

Clone the repository including submodules (this contains an adapted version of libigl that is necessary) with 

```
git clone --recursive https://github.com/FloorVerhoeven/SketchMeshVR.git
```

For the VR version, a ready-made .exe is included in bin_SketchMeshVR/bin. Simply double click it to run. If you prefer to build it yourself, follow the instructions below.


Create a build folder in the "SketchMesh" and/or "VR-SketchMesh" folder.
Use CMake to build the makefiles for your system. You can either do this via the CMake GUI (press configure and then generate; use the "Unix Makefiles" generator on macOS or the "Visual Studio 15 2017 Win64" one on Windows) or via the command line by navigating into the build folder and running:

```
cmake -DCMAKE_BUILD_TYPE=Release ../
```

Then for macOS, use the terminal to navigate to inside the build folder run the following to compile and run the executable (replace ./SketchMesh_bin with ./VR-SketchMesh_bin for the VR version:
```
make && ./SketchMesh_bin
```
This will launch the application.

## Usage of VR-SketchMesh
At all times, the control curves are overlayed over the mesh. Sharp curves are indicated in red, and smooth curves are displayed in blue.
All editing actions are performed by pressing and holding down the **TRIGGER** button (at your index finger) of the right controller.

The default action is the **DRAW** mode, which allows you to create a mesh.
When you want to start a new mesh, simply start drawing over the current one (in DRAW mode), it will be removed as soon as you start. 

Laser usage: ADD, REMOVE, CUT, EXTRUDE (first part), TYPE CHANGE
Tip-of-finger usage: DRAW, PULL, EXTRUDE (second part), SMOOTH

All of the modes below can be selected from the menu, which can be brought up by pressing the **A** button. Modes can be selected from the menu by aiming at them and then pressing either the **B** button or the **TRIGGER** button.
* **DRAW**: Press and hold down the **TRIGGER** button and draw the shape. You don't need to fully close the shape, as this will be done for you. The outline is not allowed to intersect itself. The sofware will interpret the outline as if it was drawn in 2D on a whiteboard that is facing your viewing direction.
* **ADD**: Press and hold down the **TRIGGER** button and either draw the new curve entirely on the mesh surface, or start and end outside of the mesh. When the stroke starts and ends outside of the mesh, the stroke will wrap around over the backside of the mesh (like a cut stroke).
* **REMOVE**: Point at the curve you want to remove and press the **TRIGGER** button. The selected stroke will show up in black. Point at the same stroke and press again to confirm its removal. 
* **PULL**: Hold your hand close to the curve you want to deform and press and hold down the **TRIGGER** button. Pull the curve to the desired shape by moving your hand. When your hand is too far away from any of the control curves, nothing will happen. NOTE: large deformations will result in the entire mesh being translated. To avoid this, break up large deformations into multiple smaller ones and release the **TRIGGER** button in between.
* **CUT**: Aim the laser outside of the mesh, then press and hold down the **TRIGGER** button and draw a cut stroke over the mesh as if you are cutting with a laser sword. Make sure to start and end the cut stroke outside of the mesh. After drawing the cutting stroke, point at the mesh part that you want to remove and press the **TRIGGER** button to remove it.
* **EXTRUDE**: First draw the extrusion base entirely on the mesh surface while holding down the **TRIGGER** button (it is possible to draw over existing curves). When this is successful, it will show up in blue. Then, using the position of your hand, press and hold down the **TRIGGER** button and draw the extrusion silhouette. This silhouette will be attached to the extrusion base. 
* **TYPE CHANGE**: Point at the control curve whose type you want to change from sharp to smooth or vice versa and press the **TRIGGER** button.
* **SAVE**: Selecting the save mode will open a pop-up screen on the monitor of the PC that the Oculus Rift is attached to. This pop-up asks for the location where you want to save the mesh. Use the keyboard and/or mouse to give input.
* **LOAD**: Selecting the load mode will open a pop-up screen on the monitor of the PC that the Oculus Rift is attached to. This pop-up asks for the file that you want to load a mesh from. Use the keyboard and/or mouse to give input.

In addition the the modes that can be selected in the menu, the user can also smooth the mesh or individual curves. In order to smooth the complete mesh, you press the **B** button while holding your hand away from the mesh. In order to smooth an individual control curve, hold your hand close to the curve you want to smooth, and then while holding down the **B** button rub over the curve in order to smooth it.
Finally the user can rotate the mesh in space by moving the right **thumbstick**.

## Usage of non-VR SketchMesh
At all times, the control curves are overlayed over the mesh.

* For drawing a initial mesh, press 'D' and draw the shape (while pressing the mouse down). Slower drawing results in more internal triangles but also a slower processing.
* For adding additional control curves on the existing mesh surface, press 'A' and draw the stroke (while pressing the mouse down). For strokes that start outside of the mesh but end inside it, only the points that were drawn on the mesh will be considered. Strokes that start and end off the mesh are curently not supported
* For pulling (deforming) the shape, press 'P' and click close to one of the vertices of the stroke that you want to pull on, and while holding the mouse button down, pull to the desired position. The menu on the left hand side allows for switching between smooth and sharp deformation. If you click too far away from one of the stroke vertices, the mouse movement will be interpreted as a navigation movement.
* For removing additional control curves, press 'R' and click close to one of the vertices of the stroke that you want to remove. Its vertices will turn black and you need to click one of the vertices of the same curve again to confirm the removal.
* For cutting the mesh, press 'C' and draw a cutting stroke that starts and ends outside of the mesh. After drawing the cutting stroke, click in order to perform the cut. The mesh part on the left-hand side of the drawn stroke will be removed (when looking in the direction of drawing of the stroke).
* For making extrusions, press 'E' and draw an extrusion base stroke on the mesh. This stroke needs to fully contain at least one mesh vertex. Then press 'N' to turn the mesh to get a side view in which you can draw a silhouette stroke. Press 'E' again and draw the silhouette stroke outside of the mesh.
* For navigating, press 'N' and drag.
