
from scipy.io import FortranFile

infile = FortranFile("./output/procx0_procy0/mean_intensity_0000.bin", "r")

U = infile.read_record(dtype='d')

mz=16
ny=1
nx=1
nw=1
U = U.reshape(mz,ny,nx,nw)

print(U)
