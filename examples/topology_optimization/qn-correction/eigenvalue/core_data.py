import os
import argparse
import re
from shutil import copy

p = argparse.ArgumentParser()
p.add_argument('result_folder', type=str)
p.add_argument('paropt_num', type=int)
p.add_argument('start', type=int)
p.add_argument('end', type=int)
p.add_argument('problem', type=str)
args = p.parse_args()

# Get all folders with name number-optimizer
dirs = os.listdir(args.result_folder)
r = re.compile(r"\d+-.+")
dirs = list(filter(r.match, dirs))

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
    num, omz = d.split('-')

    try:
        copy(os.path.join(args.result_folder, d, 'output_refine0.pkl'),
                os.path.join(out_folder, d))
    except:
        print("cannot copy {:s}, file might not exist.".format(os.path.join(args.result_folder, d, 'output_refine0.pkl')))

    try:
        copy(os.path.join(args.result_folder, d, 'fail.f5'),
             os.path.join(out_folder, d))
    except:
        pass

    if omz == 'paropt' or omz == 'paroptqn':
        try:
            copy(os.path.join(args.result_folder, d, 'tr_output_file0.dat'),
                 os.path.join(out_folder, d))
        except:
            pass
        try:
            copy(os.path.join(args.result_folder, d, 'output{:d}.f5'.format(args.paropt_num)),
                 os.path.join(out_folder, d))
        except:
            pass

    elif omz == 'mma':
        try:
            copy(os.path.join(args.result_folder, d, 'mma_output_file0.dat'),
                 os.path.join(out_folder, d))
        except:
            pass

    else:
        try:
            copy(os.path.join(args.result_folder, d, omz+'_output_file0.dat'),
                 os.path.join(out_folder, d))
        except:
            pass
        try:
            copy(os.path.join(args.result_folder, d, 'output_refine0.f5'),
                    os.path.join(out_folder, d))
        except:
            pass

# Copy over stdouts
for d in dirs:
    num, omz = d.split('-')

    # Compute the number of stdout
    offset = {'paropt':1, 'paroptqn':2, 'snopt':3, 'ipopt':4, 'mma':5}
    stdout_num = (int(num) - args.start)*len(offset) + offset[omz]
    stdout_name = '{:s}-n{:d}-{:d}.out-{:d}'.format(args.problem,
        args.start, args.end, stdout_num)

    try:
        copy(os.path.join(args.result_folder, stdout_name),
                os.path.join(out_folder, d))
    except:
        print("Cannot copy over stdout file: {:s}".format(stdout_name))
