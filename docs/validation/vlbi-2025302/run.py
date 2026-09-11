from pathlib import Path
import json, subprocess, sys
p=Path(__file__).resolve().parent
for name in ['audit.py','heldout.py','check_delay.py']:
    subprocess.run([sys.executable,str(p/name)],check=True,stdout=subprocess.DEVNULL)
r=json.loads((p/'report.json').read_text())
assert r['quadratic_removed_360_rows']['CIO_vs_yi_PNM_GAST']['peak_to_peak_deg'] < .001
assert r['quadratic_removed_360_rows']['ERFA_aberration_vs_yi_vlbi']['peak_to_peak_deg'] < .02
assert r['interpolation_max_deg'] < .3
assert r['paired_13_rows']['nominal_troposphere']['after_subtraction']['rms_deg'] < 2.
assert r['held_out_12_rows_fixed_prediction']['rms_deg'] < 2.
assert r['residual_delay_13_rows']['after_rms_sample'] < .002
print('Independent audit passed; report.json and phase_comparison.png updated.')
