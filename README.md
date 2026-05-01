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
```setup mrb v6_09_11```<br/>
```setup dunesw v09_78_06d00 -q e26:prof```<br/>

Then setup the module<br/>

```mrb newDev```<br/>
Run the command it tells you to run<br/>
```cd srcs```<br/>

Pull the code from GitHub<br/>
```git clone -b v09_78_upgrade https://github.com/mophoven/FD_Ereco.git```<br/>


Now we begin the build<br/>
```cd ..```<br/>
```mrb uc```<br/>
```cd ${MRB_BUILDDIR}```<br/>
```mrb z```<br/>
```mrbsetenv```<br/>
```mrb b```<br/>

Build should be successful, then to test, run on the file pre-loaded into the .fcl file<br/>

```lar -c EnergyAnalysis.fcl```<br/>

This should produce a file named FD_Ereco_CAF.root under: ```/exp/dune/app/users/<your_username>/EReco/srcs/myntuples/myntuples/MyEnergyAnalysis```<br/>

**Each subsequent login**<br/>

Setup environment inside SL7 container<br/>
```unsetup mrb```<br/>
```setup mrb v6_09_11```<br/>
```setup dunesw v09_78_06d00 -q e26:prof```<br/>
```source /exp/dune/app/users/<your_username>/EReco/localProducts_larsoft_v09_22_02_debug_e19/setup```<br/>


