# Ocotillo-py

Python version of the ocapcities used in the main code.

## Usage

Run your calculation using the `solve_opacity.py` script. Here is the usage:

```
python solve_opacity.py -h
```
```
usage: solve_opacity.py [-h] [-d DENSITY] [-t TEMPERATURE] [-n NLAMBDA] [--lambdaMin LAMBDAMIN] [--lambdaMax LAMBDAMAX] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -d DENSITY, --density DENSITY
                        Gas density in cgs units (g cm^-3). Default: 1e-9 g cm^-3
  -t TEMPERATURE, --temperature TEMPERATURE
                        Gas temperature in Kelvin
  -n NLAMBDA, --nlambda NLAMBDA
                        Integer number of wavelengths to be used. Default: nlambda = 2
  --lambdaMin LAMBDAMIN
                        Minimum wavelength in Angstroms. Default: lambdaMin = 3,000 A
  --lambdaMax LAMBDAMAX
                        Maximum wavelength in Angstroms. Default: lambdaMax = 10,000 A
  -v, --verbose         If set, print verbose output
```

## Example output

```
python solve_opacity.py
Calculating for 2 wavelength bins over the range [3000, 10000] Angstroms.

Gas thermal properties:
Temperature: 1000.00 K
Density: 1.00e-09 g cm^-3

Continuous absorption coefficients by wavelength:
Wavelength (A)     Absorp. Coeff. (1/cm)
       3000.00     9.582125765212365e-36
      10000.00     1.524825818772120e-35
          Grey     1.164118178941964e-06
```