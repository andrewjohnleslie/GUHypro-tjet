# GUHypro



What is HYPRO?
==============

HYPRO is a software that is used for Propulsion System modelling. 
HYPRO can be used to model any kind of propulsion system, thanks to its
modular structure.
It is designed to be fast in order to be easily included within
optimization loops.
Full details of HyPro with the theoretical background and the validations
are available in the [PhD thesis of Alessandro Mogavero](http://oleg.lib.strath.ac.uk/R/?func=dbin-jump-full&object_id=26892)

If you use HyPro in any work please cite the following paper as reference:

Mogavero, A., Taylor, I., and Brown, R. E., *“Hybrid Propulsion Parametric and Modular Model: a novel engine analysis tool conceived for design optimization”*,
AIAA Aviation - 19th AIAA International Space Planes and Hypersonic Systems and Technologies Conference.,
American Institute of Aeronautics and Astronautics, 06 2014, doi:[10.2514/6.2014-2787](https://pure.strath.ac.uk/portal/files/34897090/Mogavero_A_et_al_Hybrid_propulsion_parametric_and_modular_model_a_novel_engine_analysis_tool_conceived_for_design_optimization_Jun_2014.pdf).

LICENCE
=======

HyPro is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

HyPro is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with HyPro.  If not, see <http://www.gnu.org/licenses/>.

Copyright 2016 Alessandro Mogavero.

=============

GUHyPro

The above exists as unedited upon the completion of Dr Mogavero's PhD.
Since then, the project has undergone a number of developments, which have led to the creation of this separate repository, created in part as a collaboration between The Universities of Glasgow and Strathclyde.

For more information, please contact:

Mark De Luca
Room 734
James Weir South Building
University of Glasgow
Email: mark.de-luca.1@research.gla.ac.uk


=============


Prerequisites
=============

In order to compile and use HyPro, a system requires the following:

1. Linux operating system. 
    1. Ubuntu 14.04 (Recommended)
    2. Ubuntu 15.04 (Works, but with errors in OpenFOAM installation)
2. An installation of Cantera
3. An installation of Boost
4. Submodules
3. [OPTIONAL] An installation of googletest.
4. [OPTIONAL] An installation of Matlab.
5. [OPTIONAL] A C++ IDE.
    1. [CLion](https://www.jetbrains.com/clion/)
    2. [Eclipse](https://eclipse.org/downloads/)
6. Access to the source folder for the ThirdParty directory, GPC++ stored on github.com/strath-ace-labs
7. fmt now required for some test cases. Install fmt-X.X.X into ThirdParty

## 1. Installation of Linux (steps taken from How to install OpenFOAM 2.1.1 (steps taken from [here](http://www.everydaylinuxuser.com/2014/05/install-ubuntu-1404-alongside-windows.html)) ##
If you are using a windows operating system, you must follow appropriate steps in order to load Ubuntu onto your system.
Before you begin, ensure that your system is appropriately backed up.

1. Create a bootable Ubuntu USB drive

    * This will require a USB drive which is empty, as you will lose all of the data on the drive.
    * Insert the USB drive into your computer.
    * Visit http://www.ubuntu.com/download/desktop. Make sure that you choose the 64-bit version of 14.04. Click "Download" to download the file.
    * Visit http://www.pendrivelinux.com. Download the universial USB installer.
    * Start up the Universal USB installer.
    * The first thing to do is choose your distribution of choice, in this case Ubuntu 14.04, from the dropdown list.
    * Click on the "Browse" button. Find the downloaded Ubuntu ISO.
    * Select your chosen USB drive letter and make sure that the "We will format" option is checked.
    * Click "Create" to continue.
    * You should now have a bootable USB drive.
    
2. Shrink your Windows partition.

    * Windows takes up the whole of the drive when it is first installed. In order to install Ubuntu you will need to make space for it.
    * Press the "super key" (Windows key) on your keyboard and click the magnifying glass in the top right corner. In the search box start typing "Partitions".
    * Click on the option called "Create and format partitions". This will bring up the "Disk Management" screen.
    * To shrink the drive, right click on the "OS (C:)" volume and select "Shrink volume".
    * A screen will appear showing how much you can shrink the drive by. You can of course choose to shrink the drive by less than offered but never go for any more than offered as you will break your Windows 8.1 operating system if you do.
    * Click "Shrink" to continue.
    * You should now have some unallocated space.
    
3. Install Ubuntu.

    * Reboot your computer and enter the boot menu (Usually delete or F2 etc.). 
    * Select the USB drive which you have loaded Ubuntu 14.04 onto.
    * The computer will be booted into a live Ubuntu 14.04 USB image. Select the "Install Ubuntu 14.04 LTS" icon.
    * Follow the default steps, careful to choose "Install alongside Windows".
    

    
## 2. How to install Cantera (steps taken from [here](http://www.cantera.org/docs/sphinx/html/install.html#ubuntu))

1. Install the necessary packages:

    `sudo aptitude install software-properties-common`

    `sudo apt-add-repository ppa:speth/cantera`

    `sudo aptitude update`

    `sudo aptitude install cantera-python3 cantera-dev`

## 3. An installation of boost ##

1. At the time of writing, the other boost libraries required are the math toolkit and the serialization module. All boost libraries can be installed easily on Ubuntu with:

    `sudo apt-get install libboost-all-dev`

or you can try and install the required packages manually.

## 4. Submodule ##

The GPC module used in HyPro is stored as a submodule, and can be downloaded with the following process: 

You will need to request access from the University of Strathclyde to their repository for GPC++, stored on the same github account as HyPro. Please contact who initially approved your access to the HyPro repository to request this. Once access has been granted, do the following steps:

    `git init`

    `gedit .gitmodules`

Then, edit the url to: `url = https://github.com/markjohndeluca/GUgpcpp.git`. Save and quit.

    `git submodule init`

    `git submodule update`
    
To build the ThirdParty directory:

    `cd ThirdParty/gpc`

    `mkdir lib`

    `make`

** NOTE: If you have conflicting errors where the include files are all in .../c++/4.8/ but your OS is trying to build with gpp 7.0 or another more recent version, edit the "Makefile.ini" Compiler flags to the version that you are trying to use. **


## 5. [OPTIONAL] How to install googletest ##

1. Clone the google test repository from https://github.com/google/googletest.

2. Build the google test repository with:

    `cmake .`

    `make`

    `sudo make install`

3. The google testing framework should be installed on your machine.

## 6. [OPTIONAL] How to install Matlab. ##

No information for now.

## 7. [OPTIONAL] How to install an IDE of your choice. ##

### 1. CLion ###

1. Download [Clion](https://www.jetbrains.com/clion/download/#section=linux-version) from the jetbrains website.

2. Extract CLion to /opt/ for global use.

    `cd /opt/ && sudo tar -zxvf ~/Downloads/clion-*.tar.gz`
    
3. Add a symbolic link to the local bin.

    `sudo ln -s /opt/[clion-*]/bin/clion.sh /usr/local/bin/clion`
    
    * NB, you will have to find the exact path [clion-*] yourself. It will be along the lines of [clion-2016.1.2] etc.
    
4. Clion should run as follows:

    `clion`
    
### 2. Eclipse (steps taken from [here](http://ubuntuhandbook.org/index.php/2014/06/install-latest-eclipse-ubuntu-14-04/)) ###

1. Install Java.

    If you don’t have Java installed on your system. Click the link below to bring up Ubuntu Software Center and click install OpenJDK Java 7.

    [Install OpenJDK Java 7](apt://openjdk-7-jre)

2. Download Eclipse from the oracle website.

    Install the Eclipse C++ IDE from the Eclipse [website](http://www.eclipse.org/downloads/)

    You may check out your OS Type 32-bit or 64-bit by going to *System Settings -> Details -> Overview*

3. Extract Eclipse to /opt/ for global use.

    `cd /opt/ && sudo tar -zxvf ~/Downloads/eclipse-*.tar.gz`

4. Add a symbolic link to the local bin.

    `sudo ln -s /opt/eclipse/eclipse /usr/local/bin/eclipse`

5. Eclipse should run now as follows:

    `eclipse`
    
Building
========

1. Clone this repository.

    Open a terminal and type:

    `cd ${HOME}/git`
    cd 
    `git clone https://github.com/FASTT/HyPro.git` **
     
 UPDATED_METHOD**
    Create a personal access token on GitHub to use as an authentication method.
    cd
    `git clone https://<personal_token>:x-oauth-basic@github.com/FASTT/HyPro.git` 
    
    Refer to the [git user's manual](https://git-scm.com/doc) for further git commands.
    
2. Switch to the HyPro directory.

    `cd ${HOME}/git/HyPro`
    
** NOTE: Ensure you have followed the steps to download and build the required ThirdParty submodules before continuing with the build for HyPro. **

3. Create a build folder.

    `mkdir build`
    
    `cd build`
    
4. Run cmake and make

    `cmake ..`
    
    `make`
    
5. [OPTIONAL] To build the .mexa64 executables:

    `cmake -Dmatlab=ON ..`
     
    `make `
    
    
