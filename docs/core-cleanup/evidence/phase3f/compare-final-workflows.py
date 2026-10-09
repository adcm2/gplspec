from pathlib import Path
import json, math, hashlib
import numpy as np
root=Path('/tmp/gplspec-phase3-resume')
# Fixed before the final candidate is run. These are sufficient bounds for the
# Phase 3 physical RMS tolerance, not a new componentwise potential criterion.
contract={
 'potential':'max pointwise error (or harmonic addition-theorem bound) <= (1e-8+1e-12)*minimum reference magnitude lower bound. This stronger sufficient condition implies the physical weighted RMS gate for any positive weights; a failure requires the full declared gate or explicit scientific attribution.',
 'kernel':'relative unweighted L2 representation difference <= 1e-8; kernel units preserved, not potential. Independent density contraction is validated by stage3e_sensitivity.',
 'density':'pointwise 1e-12 relative + 1e-14*1342 kg/m3 physical density scale',
 'power':'max difference of sqrt(P_l/P_0) <= 1e-8; monopole-normalized degree amplitude, including near-zero degrees, bounded by the 1e-8 physical field-relative criterion',
 'geometry':'radius data 1e-12 relative + 1e-14*10000 metres; sample coordinates and indices exact',
 'known_exception':'clean_bench_10/work/Bench10/perturbed/MatrixSolution.out: inspect against recorded 3D grid-view lifetime attribution; preserve failed comparison rather than pretending baseline parity.'}
(root/'final-workflow-comparison-contract.json').write_text(json.dumps(contract,indent=2)+'\n')

def read(p):
 return np.loadtxt(p,delimiter=';',ndmin=2)
def hashfile(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def inspect(program,path,a,b):
 assert a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all()
 text=str(path);kind='potential';bound=None;lower=None
 if 'PowerSpectrum' in text:
  assert (a>=0).all() and (b>=0).all()
  d=np.sqrt(b)-np.sqrt(a);passed=np.max(abs(d))<=1e-8;kind='monopole-normalized degree amplitude';va,vb=a,b; coords=True
 elif 'PhobosRadius' in text:
  coords=np.array_equal(a[:,:2],b[:,:2]);va,vb=a[:,2],b[:,2];d=vb-va
  passed=np.all(abs(d)<=1e-12*abs(va)+1e-10);kind='radius_m'
 else:
  if 'PHOBOS_POTENTIAL_OUTER' in text:
   coords=np.array_equal(a[:,:2],b[:,:2]);va=((a[:,2]+1j*a[:,3])*(10000/3600)**2)[None,:];vb=((b[:,2]+1j*b[:,3])*(10000/3600)**2)[None,:]
   n=va.shape[1];L=math.isqrt(n)-1;assert (L+1)**2==n
   bound=math.sqrt(n/(4*math.pi))*float(np.linalg.norm(vb-va))
   lower=float(abs(va[0,0])/math.sqrt(4*math.pi)-sum(math.sqrt((2*l+1)/(4*math.pi))*np.linalg.norm(va[0,l*l:(l+1)**2]) for l in range(1,L+1)))
   kind='potential_harmonic_bound_m2_s2'
  elif 'PREM200.NO.INTEGRAL' in text:
   coords=np.array_equal(a[:,0],b[:,0]);scale=(6371000/3600)**2/math.sqrt(4*math.pi);va=a[:,1]*scale;vb=b[:,1]*scale
  elif '/lm/' in text or (program in ['clean_bench_1','clean_bench_2'] and 'ExactSolution' not in text) or ('SensitivityKernel.out' in text):
   coords=np.array_equal(a[:,0],b[:,0]);va=a[:,1::2]+1j*a[:,2::2];vb=b[:,1::2]+1j*b[:,2::2]
   if 'Sensitivity' not in text:
    n=va.shape[1];assert round(math.sqrt(n))**2==n
    bound=math.sqrt(n/(4*math.pi))*float(np.max(np.linalg.norm(vb-va,axis=1)))
    lower=float(np.min(abs(va[:,0])/math.sqrt(4*math.pi)-math.sqrt((n-1)/(4*math.pi))*np.linalg.norm(va[:,1:],axis=1)))
    kind='potential_harmonic_bound_m2_s2'
  elif 'ExactSolution' in text or 'PREM200.NO.INTEGRAL' in text:
   coords=np.array_equal(a[:,0::2],b[:,0::2]);va,vb=a[:,1::2],b[:,1::2]
  elif 'Slice' in text:
   coords=np.array_equal(a[:,0::4],b[:,0::4]) and np.array_equal(a[:,1::4],b[:,1::4]);va=a[:,2::4]+1j*a[:,3::4];vb=b[:,2::4]+1j*b[:,3::4]
  else:
   if '/Cartesian/' in text:
    assert np.array_equal(a[:,0],b[:,0]);a=a[:,1:];b=b[:,1:]
   assert a.shape[1]%5==0
   coords=all(np.array_equal(a[:,j::5],b[:,j::5]) for j in range(3))
   va=a[:,3::5]+1j*a[:,4::5];vb=b[:,3::5]+1j*b[:,4::5]
  d=vb-va
  if 'Sensitivity' in text:
   kind='kernel_representation';passed=np.linalg.norm(d)<=1e-8*np.linalg.norm(va)
  elif 'DensitySolution' in text:
   kind='density_kg_m3';passed=np.all(abs(d)<=1e-12*abs(va)+1.342e-11)
  else:
   if bound is None:bound=float(np.max(abs(d)))
   if lower is None:lower=float(np.min(abs(va)))
   passed=lower>=0 and bound<=(1e-8+1e-12)*lower
 return {'quantity':kind,'coordinates_exact':bool(coords),'values':int(va.size),'max_absolute_difference':float(np.max(abs(d))), 'relative_l2_representation_difference':float(np.linalg.norm(d)/np.linalg.norm(va)) if np.linalg.norm(va)>0 else None,'potential_weighted_rms_upper_bound':bound,'reference_magnitude_lower_bound':lower,'sufficient_scientific_bound_pass':bool(passed and coords)}
if __name__=='__main__':
 candidate=json.loads((root/'final-workflows/manifest.json').read_text());rows=[]
 for case in candidate['runs']:
  program=case['program'];base=root/('baseline-classical-workflows' if program in ['clean_bench_9','clean_bench_10'] else 'baseline-workflows')/program;new=root/'final-workflows'/program
  assert not case['missing_expected_outputs'] and not case['extra_outputs'] and case['exit_code']==0
  for path in sorted(case['outputs']):
   pa,pb=base/path,new/path
   row={'program':program,'file':path,**inspect(program,path,read(pa),read(pb)),'baseline_sha256':hashfile(pa),'candidate_sha256':hashfile(pb)};rows.append(row)
   print(program,path,row['max_absolute_difference'],row['sufficient_scientific_bound_pass'],flush=True)
 (root/'final-workflow-comparison.json').write_text(json.dumps({'baseline':'15582edab8112dbb63e0101be5b26f5135298644','criteria':contract,'files':rows},indent=2)+'\n')
