import os
import argparse
import re
from shutil import copy
from glob import glob

p = argparse.ArgumentParser()
p.add_argument('result_folder', type=str)
p.add_argument('--single_optimizer', action='store_true')
args = p.parse_args()

# Find start, end and problem
start = int(args.result_folder.split('-')[1])
end = int(args.result_folder.split('-')[2])
out_list = glob(os.path.join(args.result_folder, '*.out*'))
problem = os.path.basename(out_list[0]).split('-')[0]

# Get all folders with name number-optimizer
dirs = os.listdir(args.result_folder)
r = re.compile(r"\d+-.+")
dirs = list(filter(r.match, dirs))
rank = {'paropt':.1, 'paroptqn':.2, 'snopt':.3, 'ipopt':.4, 'mma':.5}
dirs.sort(key=lambda s:int(s.split('-')[0]) + rank[s.split('-')[1]])

# Output folder name
out_folder = args.result_folder.strip(r'/') + '-core'

# Make directories
if not os.path.isdir(out_folder):
    os.mkdir(out_folder)

for d in dirs:
    if not os.path.isdir(os.path.join(out_folder, d)):
        os.mkdir(os.path.join(out_folder, d))

# Copy over output files
for d in dirs:
    print('populating {:s}'.format(os.path.join(out_folder, d)))

    num, omz = d.split('-')

    # Copy over pickles
    try:
        for pkl in glob(os.path.join(args.result_folder, d, 'output_refine*.pkl')):
            copy(pkl, os.path.join(out_folder, d))
    except:
        print('no output_refine*.pkl in {:s}'.format(os.path.join(args.result_folder, d)))

    # Copy over failing eigenvalue f5 file
    try:
        copy(os.path.join(args.result_folder, d, 'fail.f5'),
             os.path.join(out_folder, d))
    except:
        pass

    # Copy over outputs
    if omz == 'paropt' or omz == 'paroptqn':
        prefix = 'tr_output_file'
    elif omz == 'mma':
        prefix = 'mma_output_file'
    else:
        prefix = omz+'_output_file'

    try:
        for dat in glob(os.path.join(args.result_folder, d, prefix+'*.dat')):
            copy(dat, os.path.join(out_folder, d))
    except:
        print('no {:s}*.dat in {:s}'.format(prefix, os.path.join(args.result_folder, d)))

    # Copy over f5
    if omz == 'paropt' or omz == 'paroptqn':
        try:
            outputf5 = glob(os.path.join(args.result_folder, d, 'output*.f5'))
            outputf5 = [os.path.basename(f5) for f5 in outputf5]  # Get basenames
            outputf5.sort(key=lambda x: int(x.split('.')[0][6:]))  # Sort f5 by iteration number
            copy(os.path.join(args.result_folder, d, outputf5[-1]), os.path.join(out_folder, d))

        except:
            print('failing to copy over f5 file from {:s}'.format(os.path.join(args.result_folder, d)))

    else:
        try:
            for f5 in glob(os.path.join(args.result_folder, d, 'output_refine*.f5')):
                copy(f5, os.path.join(out_folder, d))
        except:
            print('failing to copy over f5 file from {:s}'.format(os.path.join(args.result_folder, d)))

# for d in dirs:
#    num, omz = d.split('-')

    # Copy over stdouts

    offset = {'paropt':1, 'paroptqn':2, 'snopt':3, 'ipopt':4, 'mma':5}
    if args.single_optimizer:
        stdout_num = int(num) - start + 1
    else:
        stdout_num = (int(num) - start)*len(offset) + offset[omz]
    stdout_name = '{:s}-n{:d}-{:d}.out-{:d}'.format(problem,
        start, end, stdout_num)

    try:
        copy(os.path.join(args.result_folder, stdout_name),
                os.path.join(out_folder, d))
    except:
        print("Cannot copy over stdout file: {:s}".format(stdout_name))
