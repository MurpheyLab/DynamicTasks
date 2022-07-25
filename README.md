# DynamicTasks
This repository is for interfacing with the HapticMaster (ACT-3D) currently in the NeuroImaging Lab (previously in the Biodex Lab) using Windows/VS. This version renders the Ball-in-Bowl game including visualization improvements and optional joystick control for testing without the HapticMaster. This version also includes the Nail-and-Hammer game.

![games](./ExperimentGames.png)

## Connecting to the HapticMaster
Start by switching the robot on, wait for a beep, then press the green start button.
Plug the Ethernet cable into the host computer, then set up a static IP address on the host. On Windows, go to Control Panel>Network and Internet>Network Connections. Click on "Change adapter settings", then right click on "Local Area Connection". In "Properties", set the desired IPv4 IP address and subnet mask. (You should only have to do this once unless you switch to a new Ethernet adapter.)  
IP address: 10.30.203.66  
Subnet mask: 255.0.0.0  
Default gateway: 10.30.203.1  
Preferred DNS server: 10.30.203.0  
To confirm connection: `ping 10.30.203.36`  

## Setting up a new Visual Studio project that uses the HapticMaster and OpenGL
0. Copy `glut32.dll`, `glfw.dll`, and `HapticAPI2.dll` into `C:\Windows\SysWOW64` (You should only have to do this once per computer.)  
1. Create an empty project within VS  
2. Add your main file (ex: `bowlinball.cpp`) to the Source folder  
3. In Configuration Properties:  
    a. Set target to console: `Linker>System>SubSystem>Console (/SUBSYSTEM:CONSOLE)`  
    b. Turn incremental linking off: `Linker>General>Enable Incremental Linking>No (/INCREMENTAL:NO)`  
    c. Add library files: `Linker>Input>Additional Dependencies>Edit` add `HapticAPI2.lib`, `glfw3.lib`, and `glut32.lib`  
    d. Add directory(s) containing library files to search path: `VC++ Directories>Library Directories`  
    e. Add directory(s) containing header files: `VC++ Directories>Include Directories`  

## Running the code in this repo on your computer

### File Paths
If you are using this repo on a different computer, the file paths for the tasks, setup file, and textures will need to be updated. All file paths can be updated in the `paths.hpp`

### Setting Parameters
Commonly altered parameters are located in `parameters.hpp`  

### Downloading content
If the contents of `OpenGL/stb` are missing when cloning this repo, try `git submodule update --init --recursive --remote`

### Output files
Before running the games, a safe workspace has to be defined. Two files: `flags.csv` and `workspace.csv` are expected in a WorkspaceCSVs folder. All output files from the workspace setup and the games are saved to a folder called SubjectData.

## Other useful information

### HapticMaster Workspace (approximate)
x: -0.1 to 0.19  
y: -0.3 to 0.3

### Isometric Setup Note
If the isometric setup computer requires a user name, enter the FSM... domain displayed (click on "How to log in to a different domain?" to make this visible) followed by user name dewaldlab


## Data Analysis

### Frequency Decomposition
You can plot a frequency decomposition of any of the trials using this file: `./DataAnalysis/FreqAnalysis/plot_spectrum.py`

### Controls Data Analysis
`Controls_BIB.py` and `Controls_NAH.py` parses the data and stores metrics into csv files, `controls-metrics.csv` and `controls-metric-windows.csv`, which are used for the statistical analyses in R: `controls-rm-h1.r`, `controls-rm-h2-h3.r`, and `controls-rm-h4.r`. The top of each R script has run instructions.

### Stroke Data Analysis
`Stroke_BIB.py` and ??????? parses the data and stores raw energy@resonance metrics into `stroke-BIB.csv` and percent comparison metric into `stroke-BIB-percent-loss-aggregate.csv` and `stroke-BIB-percent-loss-separate.csv`.
The R scripts (add here later!) import the csv files and perform the statistical analyses. All ball-in-bowl stroke results are collected in `Stroke_BIB_results_summary.pdf`.
