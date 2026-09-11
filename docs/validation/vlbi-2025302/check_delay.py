from pathlib import Path
import json
import numpy as np
P=Path(__file__).resolve().parent
r=json.loads((P/'report.json').read_text()); d=np.load(P/'curves.npz'); obs=d['obs']; t=obs[:,0];x=(t-t.mean())/1800
vals=np.array([[-39.49279022,-45.62033463],[-39.46389008,-45.62789154],[-39.43915176,-45.62437820],[-39.41942978,-45.61465454],[-39.44995499,-45.64489365],[-39.47838974,-45.66695404],[-39.50719452,-45.68431473],[-39.50935364,-45.67153549],[-39.51205444,-45.66343307],[-39.52375412,-45.67135620],[-39.50143051,-45.65885544],[-39.48423386,-45.66839218],[-39.43239212,-45.66749191]])
res=lambda v:v-np.polynomial.polynomial.polyval(x,np.polynomial.polynomial.polyfit(x,v,2))
meas=res(vals[:,0]-vals[:,1]); nominal=res(d['nominal_troposphere'][(t/10).astype(int)])/(360*6.6e9)*1.024e9
r['residual_delay_13_rows']={'before_rms_sample':float(np.sqrt(np.mean(meas**2))), 'after_rms_sample':float(np.sqrt(np.mean((meas-nominal)**2))),'correlation':float(np.corrcoef(meas,nominal)[0,1])}
np.savetxt(P/'paired_delay_rows.txt',np.column_stack([t,vals]),header='elapsed_s yi_res_delay_sample reference_output_res_delay_sample; user supplied rows only',fmt='%.8f')
(P/'report.json').write_text(json.dumps(r,indent=2)+'\n'); print(json.dumps(r['residual_delay_13_rows'],indent=2))
