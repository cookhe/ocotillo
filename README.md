# Ocotillo
Multi-wavelength post-processing radiative transfer code. 
## setup instructions
Fork the repository to your account using the "Fork" button at the top of the webpage: https://github.com/cookhe/ocotillo

Clone the repository to your local machine with either of the following methods

SSH is the preferred method. (see github's [guide](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) to setting up ssh keys)

```
git clone git@github.com:cookhe/ocotillo.git
```

or HTTPS:

```
git clone https://github.com/cookhe/ocotillo.git
```

Add the code's directory to your `PATH` variable in your `.bashrc` file (or equivalent: e.g. `.bash_profile`, `.profile`, etc.) by adding the following lines:

```
# ocotillo code
export OCO_HOME="$HOME/<path-to-code>/ocotillo"
export PATH="$OCO_HOME/bin:$PATH"
```
where `<path-to-code>` is where you've cloned the code in the previous step. Don't forget to source the file so these changes take effect!
```
source ~/.bashrc
```
## test ocotillo

Try running the autotest.
```
rt_autotest
```
If you've set things up correctly, this can be run from anywhere. Doing so executes the file `bin/rt_autotest` that runs a few tests in `samples/` to verify the code is functional.

Each successful autotest will print out some information about the runtime.
```
 processor column time =    1.4114999999999999E-002  seconds.
 Execution time =  0.14414828431372551       micro-seconds per wavelength point per mesh point.
 
 Wall time =    3.2466000000000002E-002  seconds.
 Execution time =  0.16577818627450980       micro-seconds per wavelength point per mesh point.
```

If the test fails, you will see additional output with the difference between the diagnostic file just produced `samples/<sample name>/output/diagnostics.txt` and a reference file `samples/<sample name>/reference.out`.

```
1,10c1,10
<      1    0.15145E+12    0.10205E+13    0.35091E-26    0.35128E-26
<      2    0.15145E+12    0.10205E+13    0.35091E-26    0.35128E-26
<      3    0.15145E+12    0.10205E+13    0.35091E-26    0.35128E-26
<      4    0.15145E+12    0.10205E+13    0.35091E-26    0.35128E-26
<      5    0.94922E+14    0.22690E+15    0.11264E-11    0.39839E-12
<      6    0.94922E+14    0.22690E+15    0.11264E-11    0.39839E-12
<      7    0.15145E+12    0.10205E+13    0.35091E-26    0.35128E-26
<      8    0.15145E+12    0.10205E+13    0.35091E-26    0.35128E-26
<      9    0.15145E+12    0.10205E+13    0.35091E-26    0.35128E-26
<     10    0.15145E+12    0.10205E+13    0.35091E-26    0.35128E-26
---
>      1    0.32026E+13    0.13668E+14    0.35091E-26    0.35128E-26
>      2    0.32026E+13    0.13668E+14    0.35091E-26    0.35128E-26
>      3    0.32026E+13    0.13668E+14    0.35091E-26    0.35128E-26
>      4    0.32026E+13    0.13668E+14    0.35091E-26    0.35128E-26
>      5    0.94812E+14    0.22106E+15    0.51488E-13    0.27311E-13
>      6    0.94812E+14    0.22106E+15    0.51488E-13    0.27311E-13
>      7    0.32026E+13    0.13668E+14    0.35091E-26    0.35128E-26
>      8    0.32026E+13    0.13668E+14    0.35091E-26    0.35128E-26
>      9    0.32026E+13    0.13668E+14    0.35091E-26    0.35128E-26
>     10    0.32026E+13    0.13668E+14    0.35091E-26    0.35128E-26
```

## run on your data
Create a working directory (this can be anywhere) and link the source files.
```
mkdir rt-nest && cd rt-nest
rt_setup
```
Next, you'll need to provide information about your data set in two files `input.in` and `resolution.in`

Contents of `input.in`:
```
!  -*-f90-*-  (for Emacs)    vim:set filetype=fortran:  (for vim)
&input
  z0 = -5d15
  z1 =  5d15
  sigma_grey=0.4
  lgrey=T
  lread_athena=F
/
&temperature_input
  isoTemp = 1. !64857.04298198191
  sigma = 1e30 !2.54d14
  midTemp = 1. !61341
  floorTemp = 2000.
/
&density_input
  rho0=1d-9
  rho_floor=1d-24
  H=1.496d14
/
&gas_state_input
  fully_ionized_T = 2.d4
/
```

Contents of `resolution.in`:
```
!  -*-f90-*-  (for Emacs)    vim:set filetype=fortran:  (for vim)
integer, parameter :: nw=1   ! Number of wavelengths
integer, parameter :: nx=1,ny=1,nz=10  ! Number of grid points in computational domain
integer, parameter :: nprocx=1,nprocy=1,nprocz=1
!
integer, parameter :: nyloc=ny/nprocy,nxloc=nx/nprocx  ! Number of grid points in computational domain   
```



