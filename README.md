# FD_Ereco
**Instructions for running CAFMaker on DUNE GPVM machines**
Adapted from https://github.com/weishi10141993/myntuples

Setup your dune environment inside an SL7 container

Then cd into your app directory and make a new directory for the build

```cd /exp/dune/app/users/<your_username>```<br/>
```mkDir EReco```<br/>
```cd EReco```<br/>

Then setup mrb <br/>

```unsetup mrb```<br/>
```setup mrb v4_04_06```<br/>
```setup dunetpc v09_22_02 -q e19:debug```<br/>

Then setup the module<br/>

```mrb newDev```<br/>
Run the command it tells you to run<br/>
```cd srcs```<br/>

Pull the code from GitHub<br/>
```git clone https://github.com/mophoven/FD_Ereco.git```<br/>

Remove Plotting and information scripts inside FD_Ereco<br/>
```cd FD_EReco```<br/>
```rm -r Event\ Information\ Scripts```<br/>
```rm -r Particle\ Tracking```<br/>
```rm -r Plotting\ Macros```<br/>

Now we begin the build<br/>
```cd ..```<br/>
```mrb uc```<br/>
```cd ${MRB_BUILDDIR}```<br/>
```mrb z```<br/>
```mrb setenv```<br/>
```mrb b```<br/>

Build should be successful, then to test, run on the two files pre-loaded into the .fcl file<br/>

```lar -c EnergyAnalysis.fcl -n -2```<br/>

This should produce a file named FD_Ereco_CAF.root under: ```/exp/dune/app/users/<your_username>/EReco/srcs/myntuples/myntuples/MyEnergyAnalysis```<br/>

**Each subsequent login**<br/>

Setup environment inside SL7 container<br/>
```unsetup mrb```<br/>
```setup mrb v4_04_06```<br/>
```setup dunetpc v09_22_02 -q e19:debug```<br/>
```source /exp/dune/app/users/<your_username>/EReco/localProducts_larsoft_v09_22_02_debug_e19/setup```<br/>


