flexcalc
========

Andrew C.R. Martin, abYinformatics, November 2025
-------------------------------------------------

`flexcalc` calculates a flxibility score from an MD trajectory
supplied in a simple format:

```
      >header
      x y z
      x y z
      ...
      >header
      x y z
      x y z
      ...
      (etc)
```

The program minimizes memory usage for extremely large trajectories by
making multiple passes through the file.

The code proceeds as follows:

1. Read the file through to obtain the number of frames
2. Read a second time to calculate the average position for each atom
3. Read a third time to find the frame closest to the average
4. Read a fourth time to calculate the RMSD of each frame from the
   frame closest to the average
5. Calculate and display the average of the RMSDs

### Usage

```
   ./flexcalc trajectory-file
```

### Compiling

Assuming you have `BiopLib` installed in the standard directories (`$HOME/lib` and `$HOME/include`), you simply type:

```
   make
```

This will create the `flexcalc` executable that you can run from the
current directory (`./flexcalc`) or you can move to somewhere in your
path (such as `$HOME/bin`) and just run as `flexcalc`.

