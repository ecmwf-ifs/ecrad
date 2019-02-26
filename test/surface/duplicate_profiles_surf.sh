#!/bin/bash -e
# Duplicate profiles with building height varying from 0 to 50 metres
# You need to have the nco netcdf tools in your PATH

INPUT=$1
OUTPUT=$2

NHEIGHT=15
HEIGHT='0,1,2,3,5,7,10,15,20,25,30,35,40,45,50'

# Check for existence of NCO commands
command -v ncks >/dev/null 2>&1 || { \
 echo "###########################################################" ; \
 echo "### Error: NCO commands (ncks etc) needed but not found ###" ; \
 echo "###########################################################" ; \
 exit 1; }

T_SURF=304.2
#T_SURF=299.2
T_STREET=304.2

T_AIR=294.2
T_AIR=299.2
T_AIR=304.2
 
ncks -O --mk_rec_dmn column $INPUT tmp0.nc
ncrcat -O -n $NHEIGHT,1,0 tmp0.nc tmp_$OUTPUT
ncap2 -O \
    -s "canopy_depth(:,2)={$HEIGHT}" \
    -s "building_scale(:,2)=15" \
    -s "temperature_hl(:,137)=$T_AIR" \
    -s "skin_temperature(:,3)=$T_SURF" \
    -s "canopy_temperature(:,2)=$T_AIR" \
    -s "skin_temperature(:,4)=$T_SURF" \
    -s "skin_temperature(:,2)=$T_STREET" \
    -s "tile_fraction(:,0)=0.0" \
    -s "tile_fraction(:,1)=0.0" \
    -s "tile_fraction(:,2)=1.0" \
    tmp_$OUTPUT $OUTPUT
rm tmp0.nc tmp_$OUTPUT

#    -s "temperature_hl(:,137)=304.2" \   # Lowest air temperature
#    -s "skin_temperature(:,3)=304.2" \   # Roof temperature
#    -s "canopy_temperature(:,2)=304.2" \ # Canopy air temperature
#    -s "skin_temperature(:,4)=304.2" \   # Wall temperature
#    -s "skin_temperature(:,2)=304.2" \   # Street temperature

