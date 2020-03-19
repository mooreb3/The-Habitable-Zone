;Benjamin Moore
;17327505
;Physics and Astrophysics
;Assignment 1 - JS IDL Lab
;The Habitable Zone
;A region arpund stars where liquid water may be found


;****Part1_computing the extension of the habitable zone****


my_data=read_ascii('/cphys/ugrad/2018-19/SF/MOOREB3/Downloads/Week2_data_for_week2/data_for_week2/my_data_file.dat') ;Importing data for usage
Teff=reform(my_data.field1[1,*])  ;extracting data for T effective(K), Mass(M_0), and log(g) 
log_g=reform(my_data.field1[2,*]) ;for use in the equation found in
mass=reform(my_data.field1[0,*])  ;Habitable planets around the star Gliese 581?

l_in_sun=0.72                     ;Constants for figuring out the habitable zones
l_out_sun=1.77
G=6.67408* 1d-8
solarlum=3.826e33
a_in=2.7619e-5
b_in=3.8095e-9
a_out=1.3786e-4
b_out=1.4286e-9
Tstar=(Teff-5700)

l_in_sun=[0.72,0.95,0.76,0.51]    ;arrays used to loop through given empirical equation to give multiple limits on habitable zones, depending on varying conditions
l_out_sun=[1.77,1.67,1.95,2.40]
line=['-',':','--','-.']
name=['V-M Criterion', '0% cloud','50% cloud','100% cloud']
xpos=[0.3,0.3,0.3,0.3,0.5,0.5,0.5,0.5]
ypos=[1.0,0.8,0.6,0.4,0.8,0.6,0.4,0.2]

R=sqrt((G)*(mass*1.989e33)/(10^(log_g)))  ;radius in cm, canceling all terms except for cm^2 in bracket
Rstar=R/(6.68459* 1d10)                   ;Converting from cm to AU
lum=av_luminosity_bol(Rstar,Teff)         ;luminosity in erg/s, convert to J/s

 
;a for loop to plot the limits of the habitable zones under varying conditions of cloud cover
for i=0,3 do begin     
  l_in= (l_in_sun[i] - (a_in*Tstar) - (b_in*(Tstar^2)))*(lum/solarlum)^(1./2.)
  l_out = (l_out_sun[i] - a_out*Tstar - b_out*(Tstar^2))*(lum/solarlum)^(1./2.)
  plt=plot(l_in,mass,/overplot,linestyle=line[i], name=name[i])
  plt2=plot(l_out,mass,/overplot,linestyle=line[i], name= 'Outer Limit')
  leg=legend(target=[plt,plt2],position=[xpos[i],ypos[i], xpos[i+4],ypos[i+4]])
endfor

;setting logarithmic axis, an appropriate range, and titles
plt.xlog=1
plt.ylog=1
plt.xrange=[0.0,10]
plt.yrange=[0.07,1]
plt.title="Habitable Zone"
plt.xtitle="Orbit Distance (AU)"
plt.ytitle="Mass (solar mass)"

;****Task 2-exoplanet data****
;Exporting data from csv file, and allowing for the first 20 rows to be a header
my_data1=read_csv('/cphys/ugrad/2018-19/SF/MOOREB3/Downloads/Week2_data_for_week2/data_for_week2/exoplanets_nasa_archive_2017Oct05.csv', n_table_header=20)
mass_exoplanets=((reform(my_data1.field08[*])*1.898e30)/5.972e27) ;planet mass exported in jupiter masses, converted to earth mass
R_exo=reform(my_data1.field06[*]) ;Semi-major axis exported

star_mass=reform(my_data1.field12[*])   ;importing the stellar masses
M_accept=where(mass_exoplanets gt 0.0 and mass_exoplanets lt 5.0)   ;setting acceptance criterion for non-super earths
superearths=where(mass_exoplanets gt 5) ;criterion for super-earths, so as not to plot some values twice

star_accept=star_mass[M_accept] ;applying the criteria to the data
R_accept=R_exo[M_accept]
plt3=scatterplot(R_exo[superearths],star_mass[superearths],name='Super-Earths') ;plotting the data that passes super-earth criteria
plt4=scatterplot(R_accept,star_accept, symbol='s',sym_color='green',/sym_filled,/overplot, name='Non-super-Earth systems', title='Super-Earth Systems') ;plotting data that passes non-super-earth criterion

for i=0,3 do begin  ;for loop used for overplotting previous graph of habitable zones
  l_in= (l_in_sun[i] - (a_in*Tstar) - (b_in*(Tstar^2)))*(lum/solarlum)^(1./2.)
  l_out = (l_out_sun[i] - a_out*Tstar - b_out*(Tstar^2))*(lum/solarlum)^(1./2.)
  plt=plot(l_in,mass,/overplot,linestyle=line[i], name=name[i])
  plt2=plot(l_out,mass,/overplot,linestyle=line[i], name= 'Outer Limit')
endfor        
;setting axis titles, ranges and scales
plt4.xtitle="Orbit Distance, [AU]"
plt4.ytitle="Mass, [solar mass]"
plt3.xrange=[0.01,10]
plt3.yrange=[0.07,1]
plt4.xlog=1
plt4.ylog=1
lege=legend(target=[plt3,plt4],position=[1,0.4,1.3,0.2])  ;scatterplot legends

;****Task 3-Keplers 3rd-Period vs Orbit
period=reform(my_data1.field05[*])  ;extracting period data
accept=where(period ne 0 and R_exo ne 0)  ;setting acceptance criteria
period_accept=period[accept]  ;passing datas through acceptance criteria
R_new=R_exo[accept]

;Scatterplot of log to base 10 of the semi-major and period
newplot=scatterplot(alog10(R_new),alog10(period_accept), symbol='*', sym_color="red", title= 'log(semi major) Vs. log(period)',name='Exoplanets')
result=linfit(alog10(R_new),alog10(period_accept)) ;getting the intercept and slope of the linear fit of the data
newplot1=plot('1.4742432*(x)+2.5588323',/overplot, name='Linear fit') ;plotting a line through the same slope and intercept
newplot.xtitle="Semi_major, [AU]" 
newplot.ytitle="period, [days]"
newplot.xrange=[-3,3]
;legen=legend(target=[newplot,result],position=[1.,0.4,1.3,0.2])

print, result

;Task 4
;g/cm^3- cgs units
mass_exoplanets_cgs=(reform(my_data1.field08[*]))*1.898e30  ;mass from jupiter masses to grams (CGS)
radii_exo_cgs=reform(my_data1.field10[*])*6.9911e9          ;radii from jupiter radii to cm

criterion=where(mass_exoplanets_cgs ne 0. and radii_exo_cgs ne 0.)    ;acceptance criteria

finalplt=scatterplot(mass_exoplanets_cgs[criterion], radii_exo_cgs[criterion], xtitle='Planet Mass, [g]', ytitle= 'Radius, [cm]', Title='Mass Vs. Radius')
finalplt.magnitude=alog10((3*mass_exoplanets_cgs[criterion])/(4*!dpi*(radii_exo_cgs[criterion])^3)) ;strength of markers set
finalplt.xrange=[1d26,6d31]
finalplt.yrange=[1e1,1e10]
finalplt.xlog=1
finalplt.rgb_table=6
line1=plot('smooth((x/((4./3.)*!dpi*1.))^(1./3.),1)', /overplot, linestyle='-', name='1 g cm^-3') ;overplotting contours
line2=plot('(x/((4./3.)*!dpi*2.))^(1./3.)',/overplot, linestyle='-.',name='2 g cm^-3')
line3=plot('(x/((4./3.)*!dpi*4.))^(1./3.)',/overplot, linestyle=':', name= '4 g cm^-3')
line4=plot('(x/((4./3.)*!dpi*8.))^(1./3.)',/overplot, linestyle='__', name='8 g cm^-3')
leg=legend(target=[line1,line2,line3,line4],position=[0.3,0.8,0.5,1])
bartype=colorbar(target=finalplt.magnitude,orientation = 1, position=[0.4,0.3,.45,.8], rgb_table=6,range=[0.5,50], title='Density, [g/cm^3]', tickvalues=[0.5, 5,25, 50.], tickname=['0.05', '1',  '5', '25','50'], ticklayout=1)
;bringing in colorbar that accounts for logarithmic scales



end
