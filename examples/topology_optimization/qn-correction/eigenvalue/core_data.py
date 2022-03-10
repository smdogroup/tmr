import os
import argparse
import re
from shutil import copy
from glob import glob
from os.path import join

p = argparse.ArgumentParser()
p.add_argument('result_folder', type=str)
p.add_argument('--include_stdout', action='store_true')
p.add_argument('--stdout_optimizers', nargs='*',
    default=['paropt', 'paroptqn', 'snopt', 'ipopt', 'mma'],
    help='Make sure the order is consistent with order defined in job_array.py')
args = p.parse_args()

# Strip trailing /, if any
result_folder = args.result_folder.strip(r'/')

# Find start, end and problem
start = int(result_folder.split('-')[1])
try:
    end = int(result_folder.split('-')[2])
except:
    end = start
out_list = glob(join(result_folder, '*.out*'))
problem = os.path.basename(out_list[0]).split('-')[0]

# Get all folders with name number-optimizer
dirs = os.listdir(result_folder)
r = re.compile(r"\d+-.+")
dirs = list(filter(r.match, dirs))
rank = {'paropt':.1, 'paroptqn':.2, 'snopt':.3, 'ipopt':.4, 'mma':.5}
dirs.sort(key=lambda s:int(s.split('-')[0]) + rank[s.split('-')[1]])

# Output folder name
out_folder = result_folder + '-core'

# Make directories
if not os.path.isdir(out_folder):
    os.mkdir(out_folder)
for d in dirs:
    if not os.path.isdir(join(out_folder, d)):
        os.mkdir(join(out_folder, d))

# Copy over output files
for d in dirs:
    print('populating {:s}'.format(join(out_folder, d)))

    src = join(result_folder, d)
    dest = join(out_folder, d)

    num, omz = d.split('-')

    # Copy over pickles
    pkl_list = glob(join(src, 'output_refine*.pkl'))
    if pkl_list:
        for pkl in pkl_list:
            copy(pkl, dest)
    else:
        print('[Warning] No pkl file fonud in {:s}'.format(src))

    # Copy over failing eigenvalue f5 file
    fail_f5_list = glob(join(src, 'fail.f5'))
    fail_vtk_list = glob(join(src, 'fail.vtk'))
    if fail_f5_list and not fail_vtk_list:
        os.system('f5tovtk {:s}'.format(fail_f5_list[0]))
    fail_vtk_list = glob(join(src, 'fail.vtk'))
    if fail_vtk_list:
        copy(fail_vtk_list[0], dest)

    # Copy over failing eigenvalue f5 file
    m0_f5_list = glob(join(src, 'non_design_mass.f5'))
    m0_vtk_list = glob(join(src, 'non_design_mass.vtk'))
    if m0_f5_list and not m0_vtk_list:
        os.system('f5tovtk {:s}'.format(m0_f5_list[0]))
    m0_vtk_list = glob(join(src, 'non_design_mass.vtk'))
    if m0_vtk_list:
        copy(m0_vtk_list[0], dest)

    # Copy over outputs
    if omz == 'paropt' or omz == 'paroptqn':
        prefix = 'tr_output_file'
    elif omz == 'mma':
        prefix = 'mma_output_file'
    else:
        prefix = omz+'_output_file'

    dat_list = glob(join(src, prefix+'*.dat'))
    if dat_list:
        for dat in dat_list:
            copy(dat, dest)
    else:
        print('[Warning] No dat file fonud in {:s}'.format(src))

    # Convert f5 to vtk and copy over vtk
    # Note that we only want the last f5/vtk
    f5_list= glob(join(src, 'output*.f5'))
    vtk_list = glob(join(src, 'output*.vtk'))

    # function that gets the version number
    version_number = lambda s : int(re.findall('([0-9]+)\.', s)[-1])

    if not vtk_list and f5_list:
        f5_list.sort(key=version_number)  # Sort by version number in ascent order
        f5 = f5_list[-1]  # pick the largest one
        os.system('f5tovtk {:s}'.format(f5))

    # Now we should have vtk
    vtk_list = glob(join(src, 'output*.vtk'))
    if vtk_list:
        vtk_list.sort(key=version_number)
        vtk = vtk_list[-1]
        copy(vtk, dest)
    else:
        print('[Waning] no f5 and vtk found in {:s}'.format(src))

    # Copy over stdouts
    if args.include_stdout:
        offset = {optimizer: i+1 for i, optimizer in enumerate(args.optimizers)}
        stdout_num = (int(num) - start)*len(offset) + offset[omz]
        stdout_name = '{:s}-n{:d}-{:d}.out-{:d}'.format(problem,
            start, end, stdout_num)
        try:
            opy(join(result_folder, stdout_name),
                    join(out_folder, d))
        except:
            print("Cannot copy over stdout file: {:s}".format(stdout_name))
