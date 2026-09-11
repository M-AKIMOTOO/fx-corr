from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
r=json.loads((P/'report.json').read_text()); d=np.load(P/'curves.npz'); obs=d['obs']; t=d['t']
y=np.genfromtxt(P/'yi_model.csv',delimiter=',',names=True)
el=np.column_stack([y['el1'],y['el2']]); ztd=np.array(r['input']['nominal_ztd_m'])
mapping=1.001/np.sqrt(.002001+np.sin(el)**2)
nom=360*6.6e9*(mapping[:,1]*ztd[1]-mapping[:,0]*ztd[0])/299792458.
held=np.array([[150,135.459,-88.273],[450,210.151,-80.213],[750,262.351,-64.546],[1050,261.097,-77.052],[1350,246.872,-82.267],[1650,214.542,-91.494],[1950,216.045,-59.924],[2250,159.710,-87.266],[2550,142.888,-85.775],[2850,124.506,-107.872],[3150,199.199,-71.987],[3450,304.336,-57.306]])
ot=obs[:,0]; diff=obs[:,1]-obs[:,2]; idx=(ot/10).astype(int)
ox=(ot-1800)/1800; hx=(held[:,0]-1800)/1800
poly=np.polynomial.polynomial.polyfit(ox,diff-nom[idx],2)
held_error=held[:,1]-held[:,2]-nom[(held[:,0]/10).astype(int)]-np.polynomial.polynomial.polyval(hx,poly)
train_error=diff-nom[idx]-np.polynomial.polynomial.polyval(ox,poly)
r['held_out_12_rows_fixed_prediction']={'rms_deg':float(np.sqrt(np.mean(held_error**2))),'peak_to_peak_deg':float(np.ptp(held_error)),'errors_deg':held_error.tolist()}
r['elevation_start_end_deg']=np.rad2deg(el[[0,-1]]).tolist()
r['nominal_tropo_start_end_ns']=(nom[[0,-1]]/(360*6.6)).tolist()
r['public_sources']=['https://raw.githubusercontent.com/liberfa/erfa/master/src/c2t06a.c','https://raw.githubusercontent.com/liberfa/erfa/master/src/ab.c','https://iers-conventions.obspm.fr/conventions/content/tn36.pdf (11.9; kinematics only)','https://gssc.esa.int/navipedia/index.php/Tropospheric_Delay (4,9; nominal sensitivity)']
(P/'report.json').write_text(json.dumps(r,indent=2)+'\n')
np.savetxt(P/'paired_output_rows.txt',np.vstack([obs,held]),header='elapsed_s yi_vlbi_noacel_unwrapped_deg reference_output_unwrapped_deg; user supplied rows only',fmt='%.6f')
fig,ax=plt.subplots(2,1,figsize=(10,7),sharex=True,layout='constrained')
baseline_poly=np.polynomial.polynomial.polyfit(ox,diff,2)
physical=nom+np.polynomial.polynomial.polyval((t-1800)/1800,poly-baseline_poly)
ax[0].plot(t/60,physical,label='Nominal atmosphere (ESA public equations)',color='tab:orange')
ax[0].scatter(ot/60,diff-np.polynomial.polynomial.polyval(ox,baseline_poly),label='13 supplied output pairs',s=28,color='tab:blue',zorder=3)
ax[0].scatter(held[:,0]/60,held[:,1]-held[:,2]-np.polynomial.polynomial.polyval(hx,baseline_poly),label='12 held-out output pairs',marker='x',s=35,color='tab:green',zorder=3)
ax[0].set_ylabel('Phase difference minus quadratic [deg]'); ax[0].legend(fontsize=9); ax[0].grid(alpha=.3)
ax[1].scatter(ot/60,train_error,label='13 points (quadratic nuisance fitted)',s=28)
ax[1].scatter(held[:,0]/60,held_error,label='12 held-out points',marker='x',s=35)
ax[1].axhline(0,color='black',lw=.7); ax[1].grid(alpha=.3); ax[1].legend(fontsize=9)
ax[1].set_ylabel('After atmosphere subtraction [deg]'); ax[1].set_xlabel('Minutes from 2025-10-29 08:15 UTC')
fig.suptitle('YAMAGU32–HITACH32: independent physical-model audit at 6.6 GHz\nBlack-box output comparison; no other correlator source used')
fig.savefig(P/'phase_comparison.png',dpi=160)
print(json.dumps(r['held_out_12_rows_fixed_prediction'],indent=2))
print('troposphere start/end ns',r['nominal_tropo_start_end_ns'])
