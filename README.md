# A 1-Dimensional Ice-Pelagic-Benthic transport model, (IPBM)
## About
Coupled simulation of ice, water column, and sediment biogeochemistry.

## Supported compilers:
* recent gfortran compiler (part of GCC)
* Intel Fortran Compiler version 12.1 or higher

## How to use
At first you must have compliant compiler, [Git], [CMake], and [NetCDF] Fortran library compiled with the same Fortran compiler as used for compiling IPBM. For the VisualStudio solution under Windows pre-compiled NetCDF libraries are provided.

Then:

## Linux(bash shell):
1. Download all required programs into the one folder:

  download IPBM

  `$ git clone https://github.com/limash/IPBM.git`
   
  download [FABM] from here and switch to dev-sham branch

  `$ git clone https://github.com/limash/FABM.git --branch dev-sham`
  
  download biogeochemistry model and switch to dev-sham branch
  
  `$ git clone https://github.com/limash/brom_niva_module.git --branch dev-sham`

  and download ecosystem model

  `$ git clone https://github.com/limash/ERSEM.git`

2. Add FABMDIR and NetCDF_ROOT environment variables

  For example you can add to `~/.bashrc` current lines:

  ```
  export FABMDIR='/path/to/FABM'
  export NetCDF_ROOT='/path/to/NetCDF/bin'
  ```
  
  Don't forget reload .bashrc `$ source ~/.bashrc`

3. Make a build 

  Enter IPBM folder and execute `$ bash build_release.sh`

4. Compile the code

  From build folder execute `$ make`

5. Run IPBM

  From build folder execute `$ ./IPBM`

## Windows 10, 8:

1. Download IPBM

  Right-click in Windows Explorer within the directory where you want to place the IPBM direcrory, and choose "Git Bash Here", or use PowerShell program. In the terminal that appears type:

  `$ git clone https://github.com/limash/IPBM.git`
   
  then download [FABM] from here and switch to dev-sham branch

  `$ git clone https://github.com/limash/FABM.git --branch dev-sham`
  
  download biogeochemistry model and switch to dev-sham branch
  
  `$ git clone https://github.com/limash/brom_niva_module.git --branch dev-sham`

  and download ecosystem model

  `$ git clone https://github.com/limash/ERSEM.git`

  if using other software, use these URLs.
  
2. Add IPBMDIR environment variable (only if you are going to use pre-compiled NetCDF libraries)

  * In Search, search for and then select: System (Control Panel)
  * Click the **Advanced system settings** link
  * Click **Environment Variables**. In the section **System Variables** click **New**
  * In the **New System Variable** specify the name **IPBMDIR** and the value **path:\to\IPBM**

3. Make a build

  * Start "CMake"
  * Browse the **Where is the source code** to the **path:\to\IPBM\src**
  * Browse the **Where to build the binaries** - e.g. **path:\to\IPBM\build**
  * Click the **Configure** button. Select a build system generator, if you use Intel Visual Fortran with Visual Studio integration and want to use NetCDF libraries that come with IPBM please select a 32-bit generator.
  * Now all configuration variables for the build system are listed and you can change them according to your preferences. You need set FABM_BASE variable to the directory where you have downloaded [FABM]. Then click the **Configure** button again. Specify FABM_ERSEM_BASE as well. Then click the **Configure** button again. Select **Advanced** option and specify -DFABM_NIVA_BASE to `path/to/brom_niva_module` also.
  * Click the **Configure** button until no new (red-coloured) configuration variables appear, then press **Generate** button.

4. Compile the code

  After generating the build system, you should build the software. You can do either by opening Visual Studio and choosing **Build All** (after opening **path:\to\IPBM\build\IPBM.sln**, right click on IPBM in **Solution Explorer** and select **Set as StartUp Project**) or typing **make** if using a build system based on makefiles.

5. Run IPBM

  Now you have **IPBM.exe** file in your `path:\to\IPBM\build\Debug(Release)` directory. It needs `fabm.yaml` and `KaraSea.nc` files as input data. You can find it in `..\IPBM\data` folder. In case of running IPBM under Visual Studio remember to specify the working directory (`..\IPBM\data`).

[Git]:https://git-scm.com/downloads
[FABM]:http://fabm.net
[CMake]:https://cmake.org/
[NetCDF]:http://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html
