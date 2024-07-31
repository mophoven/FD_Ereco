# FD_Ereco
**Instructions for running CAFMaker on DUNE GPVM machines**
Adapted from https://github.com/weishi10141993/myntuples

Setup your dune environment inside an SL7 container

Then cd into your app directory and make a new directory for the build

```cd /exp/dune/app/users/<your_username>```
```mkDir EReco```
```cd EReco```

Then setup mrb

```unsetup mrb```
```setup mrb v4_04_06```
```setup dunetpc v09_22_02 -q e19:debug```

Then setup the module

```mrb newDev```
Run the command it tells you to run
```cd srcs```

Pull the code from GitHub
```git clone https://github.com/mophoven/FD_Ereco.git```

Remove Plotting and information scripts inside FD_Ereco
```cd FD_EReco```
```rm -r Event\ Information\ Scripts```
```rm -r Particle\ Tracking```
```rm -r Plotting\ Macros```

Now we begin the build
```cd ..```
```mrb uc```
```cd ${MRB_BUILDDIR}```
```mrb z```
```mrb setenv```
```mrb b```

Build should be successful, then to test, run on the two files pre-loaded into the .fcl file

```lar -c EnergyAnalysis.fcl -n -2```

This should produce a file named FD_Ereco_CAF.root under "/exp/dune/app/users/<your_username>/EReco/srcs/myntuples/myntuples/MyEnergyAnalysis"

**Each subsequent login**

Setup environment inside SL7 container
```unsetup mrb```
```setup mrb v4_04_06```
```setup dunetpc v09_22_02 -q e19:debug```
```source /exp/dune/app/users/<your_username>/EReco/localProducts_larsoft_v09_22_02_debug_e19/setup```


