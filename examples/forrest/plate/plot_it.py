import os

files = []
for k in xrange(8):
    files.append('energy_norm_solution%d.plt'%(k))

    for ks in [10, 50, 100]:
        files.append('ks%d_solution%d.plt'%(ks, k))

for f in files:
    os.system('cp %s solution.plt'%(f))
    output = f.split('.')[0] + '.png'
    os.system('tec360 -mesa -b solution_export.mcr')
    os.system('mv export.png %s'%(output))
