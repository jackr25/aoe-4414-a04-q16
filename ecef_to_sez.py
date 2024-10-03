# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
#  Text explaining script usage
# Parameters:
# o_x_km = origin x coord in km
# o_y_km = origin y coord in km
# o_z_km = origin z coord in km
# x_km = position x coord in km 
# y_km = position y coord in km
# z_km = position z coord in km
# Output:
#  outputs SEZ coordindates in order s_km e_km z_km
#
# Written by Jack Rathert
# Other contributors: None
# 
# Test Case Setup : python3 ecef_to_sez.py 822.933 -4787.187 4120.262 1131.698 -4479.324 4430.228
# import Python modules
# e.g., import math # math module
import sys # argv
import math as m

# "constants"
R_E_KM = 6378.1363
E_E = 0.081819221456

# helper functions
def calc_denom(ecc,lat_rad):
    return m.sqrt(1-ecc**2 *(m.sin(lat_rad))**2)

## function description
class numpy_lite:
    def matrix_mult(self,matrix1, matrix2):
        rowsA = len(matrix1)
        colsA = len(matrix1[0])
        rowsB = len(matrix2)
        colsB = len(matrix2[0])
        if colsA != rowsB:
            return "Cannot Multiply!"
        else:
            result = [[0 for row in range(colsB)] for col in range(rowsA)]
            for i in range(rowsA):
                for j in range(colsB):
                    for k in range(colsA):
                        result[i][j] += matrix1[i][k] * matrix2[k][j]
            return result
    def matrix_add(list1,list2):
        return [x+y for x,y in zip(list1, list2)]
    def matrix_sub(list1,list2):
       return [x-y for x,y in zip(list1,list2)]
## Should have let me use numpy man
npl = numpy_lite()

def calc_ecef_to_llh(r_x_km, r_y_km, r_z_km):
  
  lat_rad = m.asin(r_z_km/m.sqrt(r_x_km**2 + r_y_km**2 + r_z_km**2))
  lon_rad = m.atan(r_y_km/r_x_km)
  r_lon_km = m.sqrt(r_x_km**2 + r_y_km**2)
  # r_z_km = r_z_km
  prev_lat_rad = float('nan')
  count = 0
  while (m.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad>10e-12)) and count<10:
   denom = calc_denom(E_E, lat_rad)
   c_e = R_E_KM/denom
   s_e = R_E_KM*(1 - E_E**2)/denom
   prev_lat_rad = lat_rad
   
   lat_rad = m.atan(( r_z_km + c_e*E_E**2 * m.sin(lat_rad) ) / r_lon_km)
   count = count + 1
  
  hae_km = r_z_km/m.sin(lat_rad) - s_e
  r_lon_km = (c_e + hae_km) * m.cos(lat_rad)
  r_z_km = (s_e+hae_km)*m.sin(lat_rad)

  print(lat_rad)
  print(lon_rad)
  print(hae_km)
  return [lat_rad, lon_rad, hae_km]


# initialize script arguments
o_x_km = float('nan')
o_y_km = float('nan')
o_z_km = float('nan')
x_km = float('nan')
y_km = float('nan')
z_km = float('nan')

# parse script arguments
if len(sys.argv)==7:
  o_x_km = float(sys.argv[1])
  o_y_km = float(sys.argv[2])
  o_z_km = float(sys.argv[3])
  x_km = float(sys.argv[4])
  y_km = float(sys.argv[5])
  z_km = float(sys.argv[6])
  ...
else:
  print(\
   'Usage: '\
   'python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km'\
  )
  exit()

# write script below this line
ecef_id = [x_km, y_km, z_km]
ecef_o = [o_x_km, o_y_km, o_z_km]

originvec = calc_ecef_to_llh(o_x_km, o_y_km, o_z_km)

idovec_ecefi =[ x-y for x,y in zip(ecef_id,ecef_o)]

idovec_ecef = [[idovec_ecefi[0]],[idovec_ecefi[1]],[idovec_ecefi[2]]]

phi = originvec[0]
theta = originvec[1]

rz_inv = [[m.sin(phi),0,-m.cos(phi)],[0,1,0],[m.cos(phi),0,m.sin(phi)]]
ry_inv = [[m.cos(theta),m.sin(theta),0], [-m.sin(theta),m.cos(theta),0],[0,0,1]]

rsezi = npl.matrix_mult(rz_inv,ry_inv)
rsez = npl.matrix_mult(rsezi,idovec_ecef)

print(rsez[0][0])
print(rsez[1][0])
print(rsez[2][0])