import sys, os
import numpy as np

np.savetxt('%iGeV_%ish_%ideg_%i%s_hist_%s.dat'%(E_prim,nshower,theta_p,pointing_angle,pointing,type_of_hist)
  ,np.transpose([mids, hist_c[0], hist_c[1], hist_f[0], hist_f[1]])
  ,newline = '\n'
  ,fmt = "%7.2f %1.6e %1.6e %1.6e %1.6e"
  ,header = (' Num_showers:%i \n E_primary (GeV): %i \n ID_prim_particle: %s \n Seeds: %i, %i \n'
              %(nshower,E_prim,prim_part,seed1,seed2) +
             ' Theta prim. part. incidence: %i deg \n Obs level (m): %i \n Atmosp model: %i'
              %(theta_p,obs_level ,atm_mod) +
             '\n Cerenk_bunch_size: %i \n Fluor_bunch_size: %i'
              %(cersize,fluorsize) +
             '\n  \n Distance to shower axis (m) | Phot_density_Cher/fluor (1/m2)'
             )
            )
print('Histogram stored into: %iGeV_%ish_%ideg_%i%s_hist_%s.dat'
       %(E_prim,nshower,theta_p,pointing_angle,pointing,type_of_hist)
      )