"""Independent yi-corr diagnostic. No other correlator implementation is used.

Inputs: previously supplied observation conditions (reconstructed, not original XML).
Reference: ERFA C2t06a/Ab; IERS TN36 Eq 11.9 kinematics (gravity excluded);
ESA Navipedia Tropospheric Delay Eqs 4,9 for a nominal sensitivity experiment.
Black-box comparison: 13 paired rows transcribed from user-provided output only.
"""
from pathlib import Path
import io, json, subprocess
import numpy as np
import erfa

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
lib=next((ROOT/'target/debug/deps').glob('liberfa-*.rlib'))
subprocess.run(['rustc','--edition=2021',str(OUT/'probe.rs'),'-L',f'dependency={ROOT}/target/debug/deps','--extern',f'erfa={lib}','-O','-o',str(OUT/'probe')],check=True)
raw=subprocess.check_output([str(OUT/'probe')],text=True)
(OUT/'yi_model.csv').write_text(raw)
y=np.genfromtxt(io.StringIO(raw),delimiter=',',names=True)
t=y['t']; mjd=60977.34375+t/86400.; tt=mjd+69.184/86400.
C=299792458.; F=6.6e9; AU=149597870700.
a=np.array([-3502544.587,3950966.235,3566381.192]); b=np.array([-3961788.974,3243597.492,3790597.692])
ra=np.deg2rad((17+33/60+2.70628/3600)*15); dec=-np.deg2rad(13+4/60+49.5482/3600)
k=erfa.s2c(ra,dec)
# C ERFA CIO-based transform, independent from yi-corr's Rust equinox path.
rot=erfa.c2t06a(2400000.5,tt,2400000.5,mjd,0.,0.)
base=np.einsum('nji,j->ni',rot,b-a)
r2=np.einsum('nji,j->ni',rot,b)
geo=-base@k/C
pvh,pvb=erfa.epv00(2400000.5,tt)
v=pvb['v']*AU/86400.; beta=v/C
proper=erfa.ab(k,beta,np.linalg.norm(pvh['p'],axis=1),np.sqrt(1.-np.sum(beta**2,axis=1)))
ab=-np.sum(base*proper,axis=1)/C
# Kinematic part ONLY of IERS (11.9). This is not a complete delay standard.
w2=np.einsum('nji,j->ni',rot,np.cross([0.,0.,7.29211514670698e-5],b))
sv=v@k/C; vb=np.sum(v*base,axis=1)/C**2
kin=(geo*(1.-np.sum(v*v,axis=1)/(2*C**2)-np.sum(v*w2,axis=1)/C**2)-vb*(1.+sv/2.))/(1.+sv+w2@k/C)
# Nominal atmosphere: independently sourced mapping and zenith model from ESA.
# Treat ellipsoid height as approximate sea-level height for sensitivity only.
h=np.array([erfa.gc2gd(1,z)[2] for z in (a,b)])
ztd=2.3*np.exp(-.116e-3*h)+.1
el=np.column_stack([y['el1'],y['el2']])
mapping=1.001/np.sqrt(.002001+np.sin(el)**2)
tropo=(mapping[:,1]*ztd[1]-mapping[:,0]*ztd[0])/C
assert np.max(np.abs(tropo-y['troposphere_s'])) < 1e-16, 'Rust atmosphere differs from independent model'

def detrend(time,values):
    x=(time-time.mean())/1800.
    return values-np.polynomial.polynomial.polyval(x,np.polynomial.polynomial.polyfit(x,values,2))
def stat(values):
    return dict(rms_deg=float(np.sqrt(np.mean(values**2))),peak_to_peak_deg=float(np.ptp(values)),start_deg=float(values[0]),end_deg=float(values[-1]))
components={
    'CIO_vs_yi_PNM_GAST':360*F*(geo-y['pnm_geo']),
    'ERFA_aberration_vs_yi_vlbi':360*F*(ab-y['pnm_minus']),
    'IERS_kinematics_vs_yi_vlbi':360*F*(kin-y['pnm_minus']),
    'VLBI_preset_correction_vs_old':-360*F*(y['pnm_minus']-y['mean_anchored']),
    'nominal_troposphere':360*F*tropo,
}
res={n:detrend(t,v) for n,v in components.items()}
report={'input':{'epoch_mjd':60977.34375,'ra_rad':ra,'dec_rad':dec,'ant1_m':a.tolist(),'ant2_m':b.tolist(),'carrier_hz':F,'tt_minus_utc_s':69.184,'dut1_xp_yp':0.,'nominal_ztd_m':ztd.tolist()},
        'quadratic_removed_360_rows':{n:stat(v) for n,v in res.items()},
        'interpolation_max_deg':float(np.max(y['interp_error_s'])*360*F)}
# Every row is directly copied from the user's printed post-correction phases.
# Unwrap full series first; the values below preserve that branch.
obs=np.array([[0,115.220,-62.499],[300,174.236,-86.795],[600,219.327,-92.949],[900,290.525,-44.924],[1200,253.251,-82.624],[1500,232.970,-85.906],[1800,210.779,-80.486],[2100,198.836,-61.988],[2400,156.151,-79.630],[2700,130.191,-96.860],[3000,163.147,-83.271],[3300,211.878,-96.980],[3590,351.965,-74.596]])
ot=obs[:,0]; idx=(ot/10).astype(int); diff=obs[:,1]-obs[:,2]
od=detrend(ot,diff)
report['paired_13_rows']={'observed_difference':stat(od)}
for n,v in components.items():
    r=detrend(ot,v[idx]); report['paired_13_rows'][n]={'after_subtraction':stat(od-r),'correlation':float(np.corrcoef(od,r)[0,1])}
report['paired_13_rows']['troposphere_free_scale']=float(od@detrend(ot,components['nominal_troposphere'][idx]) / np.sum(detrend(ot,components['nominal_troposphere'][idx])**2))
(OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
np.savez(OUT/'curves.npz',t=t,obs=obs,observed_detrended=od,**res)
print(json.dumps(report,indent=2))
