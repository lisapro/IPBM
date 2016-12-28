# BROM2
## About
A coupled benthic-pelagic model for simulation of ice, water and sediment biogeochemistry. Now it has only 1 test case for Kara sea area.

## Supported compilers:
* recent gfortran compiler (part of GCC)
* Intel Fortran Compiler version 12.1 or higher

## How to use
At first you must download [FABM] and do all prerequisites it needs (you should have compliant compiler, [Git], [CMake], and [NetCDF] Fortran library compiled with the same Fortran compiler as used for compiling BROM2 & FABM. For the VisualStudio solution under Windows pre-compiled NetCDF libraries are provided.)

Check **<FABMDIR>/src/drivers/brom/fabm_driver.h** (<FABMDIR> is the root of the FABM source tree). It must be tuned like 1d model:
```
#define _FABM_DIMENSION_COUNT_ 1
#define _FABM_DEPTH_DIMENSION_INDEX_ 1
#define _FABM_VECTORIZED_DIMENSION_INDEX_ 1
```

Then:

## Linux(bash shell):
1. Download BROM2

  `$ git clone https://github.com/limash/BROM2.git`
  
  download biogeochemistry module from the same directory as BROM2
  
  `$ git clone https://github.com/limash/brom_niva_module.git`

2. Add FABMDIR and NetCDF_ROOT environment variables

  For example you can add to `~/.bashrc` current lines:

  ```
  export FABMDIR='/path/to/FABM'
  export NetCDF_ROOT='/path/to/NetCDF/bin'
  ```
  
  Don't forget reload .bashrc `$ source ~/.bashrc`

3. Make a build 

  Enter BROM2 folder and execute `$ bash build.sh`

  or

  You can greate build folder manually: creat it, copy there files from `/data` folder, then execute command:

  `$ cmake path/to/BROM2/code -DFABM_BASE=$FABMDIR -DFABM_NIVA_BASE=path/to/brom_niva_module`

  from it.

4. Compile the code

  From build folder execute `$ make`

5. Run BROM2

  From build folder execute `$ ./brom2`

## Windows 10, 8:

1. Download BROM2

  Right-click in Windows Explorer within the directory where you want to place the BROM2 direcrory, and choose "Git Bash Here", or use PowerShell program. In the terminal that appears type:

  `$ git clone https://github.com/limash/BROM2.git`
  
  and to download biogeochemistry module type
  
  `$ git clone https://github.com/limash/brom_niva_module.git`

  if using other software, use these URLs.
  
2. Add BROMDIR environment variable

  * In Search, search for and then select: System (Control Panel)
  * Click the **Advanced system settings** link
  * Click **Environment Variables**. In the section **System Variables** click **New**
  * In the **New System Variable** specify the name **BROMDIR** and the value **path:\to\BROM2**

3. Make a build

  * Start "CMake"
  * Browse the **Where is the source code** to the **path:\to\BROM2\src**
  * Browse the **Where to build the binaries** - e.g. **path:\to\BROM2\build**
  * Click the **Configure** button. Select a build system generator, if you use Intel Visual Fortran with Visual Studio integration and want to use NetCDF libraries that come with BROM2 please select a 32-bit generator.
  * Now all configuration variables for the build system are listed and you can change them according to your preferences. You need set FABM_BASE variable to the directory where you have downloaded [FABM], also you need specify -DFABM_NIVA_BASE to `path/to/brom_niva_module`.
  * Click the **Configure** button until no new (red-coloured) configuration variables appear, then press **Generate** button.

4. Compile the code

  After generating the build system, you should build the software. You can do either by opening Visual Studio and choosing **Build All**, or typing **make** if using a build system based on makefiles.

5. Run BROM2

  Now you have **brom2.exe** file in your `path:\to\BROM2\build\Debug(Release)` directory. It needs `fabm.yaml` and `KaraSea.nc` files as input data. You can find it in `..\BROM2\data` folder. In case of running BROM2 under Visual Studio remember to specify the working directory (`..\BROM2\data`).

[Git]:https://git-scm.com/downloads
[FABM]:http://fabm.net
[CMake]:https://cmake.org/
[NetCDF]:http://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html
