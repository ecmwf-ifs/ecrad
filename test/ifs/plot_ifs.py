from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

#% This python scripts plots some of the inputs to the radiation scheme
#% in Fig. 1 and the outputs in Fig. 2
code = 'ecrad_meridian'
in_file=Dataset(code+'.nc','r')
pressure_hl=in_file.variables['pressure_hl'][:]
temperature_hl=in_file.variables['temperature_hl'][:]
lat=in_file.variables['lat'][:]
cos_solar_zenith_angle=in_file.variables['cos_solar_zenith_angle'][:]
skin_temperature=in_file.variables['skin_temperature'][:]
cloud_fraction=in_file.variables['cloud_fraction'][:]
aerosol_mmr=in_file.variables['aerosol_mmr'][:]
cases = [code+'_noaer_out.nc',
	 code+'_default_out.nc',
	 code+'_expran_out.nc',
	 code+'_tc_out.nc',
	 code+'_spartacus_out.nc',
	 code+'_spartacus_maxentr_out.nc']


f = [0]*6
k=0
for icase in cases:
    f[k] = Dataset(icase,'r')
    k=k+1



case_list = [0,1,2,3,4,5]
leg = ['McICA Exp-Exp no aerosols',
       'McICA Exp-Exp',
       'McICA Exp-Ran',
       'Tripleclouds Exp-Ran',
       'SPARTACUS Exp-Ran',
       'Classic SPARTACUS Exp-Ran']

styles = ['b','r','g','m','c','k']


p = 0.01*0.5*np.median(pressure_hl[:,0:-1]+pressure_hl[:,1::],0)

nplot=4
plt.figure(1,figsize=(10,10))
# set(gcf,'defaultlinelinewidth',1);
# set(gcf,'units','inches','paperposition',[0.5 0.5 15 30]);
# nplot = 4;
#fig, (ax0, ax1) = plt.subplots(nrows=2)

plt.subplot(nplot,1,1)
plt.plot(lat,cos_solar_zenith_angle,color='k', linewidth=2.0)
ax= plt.gca()
ax.set_xlim([-90,90])
plt.ylabel('Cos solar zenith ang.')


plt.subplot(nplot,1,2)
plt.plot(lat,skin_temperature-273.15,color='r',linewidth=2.0,label='Skin temperature')
plt.plot(lat,temperature_hl[:,-1]-273.15,color='b', linewidth=2.0,label='Air between first two model levels')
plt.gca().set_xlim([-90,90])
plt.ylabel('Temperature ($^\circ$C)');
plt.gca().set_ylim([-50,60]);
plt.gca().legend(loc='best',fontsize='xx-small')

plt.subplot(nplot,1,3)
CS=plt.contourf(lat,p,cloud_fraction.T,10, vmin=0.05, vmax=0.95,cmap="viridis")
plt.colorbar()
plt.contour(CS,levels=CS.levels[::4],colors='k')
plt.gca().set_ylim([0,1013])
plt.gca().invert_yaxis()
plt.gca().set_xlim([-90,90])
plt.ylabel('Pressure (hPa)')
textstyle = dict(size=10, color='white')
plt.gca().text(-70,100,'Cloud fraction',**textstyle)
# text(-90, 0, [' ' 10 '  Cloud fraction'],'verticalalignment','top')


print(aerosol_mmr.shape)
plt.subplot(nplot,1,4)
sea_salt = 1e9*np.sum(aerosol_mmr[:,0:2,-1],1);
dust = 1e9*np.sum(aerosol_mmr[:,3:5,-1],1);
organics = 1e9*np.sum(aerosol_mmr[:,6:7,-1],1);
black_carbon = 1e9*np.sum(aerosol_mmr[:,8:9,-1],1);
sulphate = 1e9*np.sum(aerosol_mmr[:,10:11,-1],1);

plt.semilogy(lat,sea_salt,'b',label='Sea salt')
plt.semilogy(lat,dust,'r',label='Dust')
plt.semilogy(lat,organics,'g',label='Organics')
plt.semilogy(lat,black_carbon,'k',label='Black carbon')
plt.semilogy(lat,sulphate,'m',label='Sulphate')
plt.gca().legend(loc='best',fontsize='xx-small')
plt.ylabel('Aerosol mass mixing ratio\n ($\mu$g/kg)')
plt.xlabel('Latitude (deg N)')
plt.gca().set_xlim([-90,90])

plt.savefig('input.png')

nplot=4
plt.figure(2,figsize=(10,10))
# set(gcf,'defaultlinelinewidth',1);
# set(gcf,'units','inches','paperposition',[0.5 0.5 15 30]);

plt.subplot(nplot,1,1)
for icase in case_list:
    print(leg[icase])
    plt.plot(lat, f[icase].variables['flux_dn_sw'][:,-1],label=leg[icase])
    plt.gca().set_xlim([-90, 90]);
    plt.ylabel('Surface SW down \n(W m$^{-2}$)');
    plt.gca().legend(loc='best',fontsize='xx-small')


plt.subplot(nplot,1,2)
for icase in case_list:
    cre =  f[icase].variables['flux_dn_sw'][:,1]-f[icase].variables['flux_dn_sw_clear'][:,1]-f[icase].variables['flux_up_sw'][:,1]+f[icase].variables['flux_up_sw_clear'][:,1]
    plt.plot(lat, cre,label=leg[icase])
    plt.gca().set_xlim([-90, 90]);
    plt.ylabel('SW cloud radiative effect \n(W m$^{-2}$)');
    plt.gca().legend(loc='best',fontsize='xx-small')


plt.subplot(nplot,1,3)
for icase in case_list:
    cre =  f[icase].variables['flux_dn_lw'][:,1]-f[icase].variables['flux_dn_lw_clear'][:,1]-f[icase].variables['flux_up_lw'][:,1]+f[icase].variables['flux_up_lw_clear'][:,1]
    plt.plot(lat, cre,label=leg[icase])
    plt.gca().set_xlim([-90, 90]);
    plt.ylabel('LW cloud radiative effect \n(W m$^{-2}$)');
    plt.gca().legend(loc='best',fontsize='xx-small')

plt.subplot(nplot,1,4)
for icase in case_list:
    plt.plot(lat, f[icase].variables['cloud_cover_sw'],label=leg[icase]);
    plt.gca().set_xlim([-90, 90]);
    plt.ylabel('Cloud cover');
    plt.gca().legend(loc='best',fontsize='xx-small')
    plt.xlabel('Latitude (deg N)');


plt.savefig('output.png')
#plt.show()
