import os
import argparse
from shutil import copy

p = argparse.ArgumentParser()
p.add_argument('result-folder', type=str)
p.add_argument('paropt-num', type=int)
args = p.parse_args()

# Get all folders with name e-number-optimizer
dirs = os.listdir(args.result_folder)
r = re.compile(r"e-\d+-.+")
dirs = list(filter(r.match, dirs))

# Output folder name
out_folder = args.result_folder + '-core'

# Make directories
if not os.path.isdir(out_folder):
    os.mkdir(out_folder)

for d in dirs:
    if not os.path.isdir(os.path.join(out_folder, d)):
        os.mkdir(os.path.join(out_folder, d))

# Copy over files
for d in dirs:
    e, num, omz = d.split('-')

    copy(os.path.join(args.result_folder, d, 'output_refine0.pkl'),
            os.path.join(out_folder, d))

    if omz == 'paropt' or omz == 'paroptqn':
        copy(os.path.join(args.result_folder, d, 'tr_output_file0.dat'),
             os.path.join(out_folder, d))
        copy(os.path.join(args.result_folder, d, 'output{:d}.f5'.format(args.paropt_num)),
             os.path.join(out_folder, d))

    else:
        copy(os.path.join(args.result_folder, d, omz+'_output_file0.dat'),
             os.path.join(out_folder, d))
        copy(os.path.join(args.result_folder, d, 'output_refine0.f5'),
                os.path.join(out_folder, d))




