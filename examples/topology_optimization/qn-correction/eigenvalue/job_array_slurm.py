from csv import DictReader
import argparse
import os
from pprint import pprint


def dic2str(dc):
    '''
    Convert a dictionary to a string with argument-type format
    example:
    dc = {'no':1, 'domain':'cantilever', 'AR': 1}
    dic2str(dc) = '--domain cantilever --AR 1'
    
    Note:
    key name starting with 1 or more underscores
    are considered comments and will not be included
    '''
    row = ['--{:s} {:s}'.format(k, v) for k, v in dc.items() if k[0] != '_']
    return ' '.join(row[1:])

def optimizer2cmd(optimizer, problem, compfreq_qn):
    if optimizer == 'paroptqn':
        if problem == 'compfreq':
            if compfreq_qn == 'comp':
                return '--optimizer paropt --qn-correction-comp'
            elif compfreq_qn == 'freq':
                return '--optimizer paropt --qn-correction-freq'
            elif compfreq_qn == 'compfreq':
                return '--optimizer paropt --qn-correction-comp --qn-correction-freq'
        else:
            return '--optimizer paropt --qn-correction'
    else:
        return '--optimizer {:s}'.format(optimizer)

def readCSV(csvfile, start, end):
    # Create csv dictionary reader
    with open(csvfile, mode='r', encoding='utf-8-sig') as f:
        reader = DictReader(f)

        # Load case dictionaries in a list
        physics = []
        for row in reader:
            no = int(row['no'])
            if no >= start and no <= end:
                physics.append(row)
    return physics

def createCaseCmds(problem, physics, optimizers, exescript, compfreq_qn):
    # Create runcase commands
    case_cmds = []
    n_exist_case = 0

    for case_dict in physics:
        for optimizer in optimizers:

            # Name of case folder
            # Example: 1-paropt
            case_folder = '{:s}-{:s}'.format(case_dict['no'], optimizer)

            # If case is not done already, create cmd
            try:
                n_mesh_refine = int(case_dict['n-mesh-refine'])
            except KeyError:
                n_mesh_refine = 1
            pkl_name = 'output_refine{:d}.pkl'.format(n_mesh_refine - 1)

            if not os.path.isfile(os.path.join(case_folder, pkl_name)):
                cmd = 'python {:s} {:s} {:s} --prefix {:s}'.format(
                    exescript,                         # eig-max.py
                    dic2str(case_dict),                # --domain cantilever --AR 1.0 ...
                    optimizer2cmd(optimizer, problem, compfreq_qn), # --optimizer paropt --qn-correction
                    case_folder                        # --prefix 1-paroptqn
                )
                case_cmds.append(cmd)
            else:
                n_exist_case += 1

    return case_cmds, n_exist_case

def writeCmdToTxt(textfile, case_cmds):
    with open(textfile, mode='w') as f:
        for row in case_cmds:
            f.write(row + '\n')
    return

if __name__ == '__main__':

    # Argument parser
    p = argparse.ArgumentParser()
    p.add_argument('--start', type=int, default=1)
    p.add_argument('--end', type=int, default=20)
    p.add_argument('--problem', type=str, default='eig',
        choices=['eig', 'comp', 'freq', 'compfreq'])
    p.add_argument('--compfreq-qn', type=str, default='comp',
        choices=['comp', 'freq', 'compfreq'])
    p.add_argument('--optimizer', type=str, nargs='*',
        default=['paropt', 'paroptqn', 'snopt', 'ipopt', 'mma', 'mma4py'],
        choices=['paropt', 'paroptqn', 'snopt', 'ipopt', 'mma', 'mma4py'])
    p.add_argument('--walltime', type=int, default=24, help='in hours')
    args = p.parse_args()

    if args.problem == 'eig':
        csv = 'eig_cases_main.csv'
        output = '_eig_cases.txt'
        pbs = '_eig_cases.sbatch'
        exescript = 'eig-max.py'
    elif args.problem == 'comp':
        csv = 'comp_cases_main.csv'
        output = '_comp_cases.txt'
        pbs = '_comp_cases.sbatch'
        exescript = 'comp-min.py'
    elif args.problem == 'freq':
        csv = 'freq_constr.csv'
        output = '_freq_cases.txt'
        pbs = '_freq_cases.sbatch'
        exescript = 'frequency.py'
    elif args.problem == 'compfreq':
        csv = 'compfreq_cases_main.csv'
        output = '_compfreq_cases.txt'
        pbs = '_compfreq_cases.sbatch'
        exescript = 'comp-min-freq-constr.py'

    # Create list of physics problems
    physics = readCSV(csv, args.start, args.end)

    # Create case commands
    case_cmds, n_exist_case = createCaseCmds(args.problem, 
        physics, args.optimizer, exescript, args.compfreq_qn)

    # Write case commands to txt
    writeCmdToTxt(output, case_cmds)

    # Print out summary
    n_planned_cases = len(args.optimizer)*len(physics)
    n_new_cases = len(case_cmds)
    print("------ Summary ------")
    print("Total number of cases needed:{:4d}".format(n_planned_cases))
    print("Number of existing cases:    {:4d}".format(n_exist_case))
    print("Number of new cases created: {:4d}".format(n_new_cases))

    if(n_new_cases + n_exist_case != n_planned_cases):
        raise RuntimeError('Case numbers don\'t match!')

    # Generate pbs
    with open(pbs, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('#SBATCH --job-name={:s}-n{:d}-{:d}\n'.format(args.problem, args.start, args.end))
        f.write('#SBATCH --account=gts-gkennedy9-coda20\n')
        f.write('#SBATCH --nodes=1 --ntasks-per-node=24\n')
        f.write('#SBATCH --mem-per-cpu=6gb\n')
        f.write('#SBATCH --time={:d}:00:00\n'.format(args.walltime))
        f.write('#SBATCH --output=%a.out\n')
        f.write('#SBATCH --array=1-{:d}\n'.format(n_new_cases))
        f.write('\n')
        f.write('cd $SLURM_SUBMIT_DIR\n')
        f.write('\n')
        f.write('#Get ith line of target file and run\n')
        f.write('export casecmd=`cat {:s} | head -$SLURM_ARRAY_TASK_ID | tail -1`\n'.format(output))
        f.write('\n')
        f.write('srun $casecmd\n')
        f.write('\n')

