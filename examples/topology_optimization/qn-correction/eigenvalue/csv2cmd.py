from csv import DictReader
import argparse
import os


def dic2str(dc):
    '''
    Convert a dictionary to a string with argument-type format
    example:
    dc = {'no':1, 'domain':'cantilever', 'AR': 1}
    dic2str(dc) = '--domain cantilever --AR 1'
    '''
    row = ['--{:s} {:s}'.format(k, v) for k, v in dc.items()]
    return ' '.join(row[1:])

def optimizer2cmd(optimizer):
    if optimizer == 'paroptqn':
        return '--optimizer paropt --qn-correction'
    else:
        return '--optimizer {:s}'.format(optimizer)

if __name__ == '__main__':

    # Argument parser
    p = argparse.ArgumentParser()
    p.add_argument('--csv', type=str, default='eig_cases.csv')
    p.add_argument('--output', type=str, default='eig_cases.txt')
    p.add_argument('--optimizer', type=str, nargs='*',
        default=['paropt', 'paroptqn', 'snopt', 'ipopt'],
        choices=['paropt', 'paroptqn', 'snopt', 'ipopt'])
    p.add_argument('--exec', type=str,
        default='~/git/tmr/examples/topology_optimization/qn-correction/eigenvalue/eig-max.py')
    p.add_argument('--walltime', type=int, default=24)
    args = p.parse_args()

    # Create csv dictionary reader
    with open(args.csv, mode='r', encoding='utf-8-sig') as f:
        reader = DictReader(f)

        # Load case dictionaries in a list
        cases = []
        for row in reader:
            cases.append(row)

    # Create runcase commands
    case_cmds = []
    n_exist_case = 0

    for case_dict in cases:
        for optimizer in args.optimizer:

            # Name of case folder
            # Example: e-1-paropt
            case_folder = 'e-{:s}-{:s}'.format(case_dict['no'], optimizer)

            # If case is not done already, create cmd
            n_mesh_refine = int(case_dict['n-mesh-refine'])
            pkl_name = 'output_refine{:d}.pkl'.format(n_mesh_refine - 1)

            if not os.path.isfile(os.path.join(case_folder, pkl_name)):
                cmd = 'python {:s} {:s} {:s} --prefix {:s}'.format(
                    args.exec,                         # eig-max.py
                    dic2str(case_dict),                # --domain cantilever --AR 1.0 ...
                    optimizer2cmd(optimizer),          # --optimizer paropt --qn-correction
                    case_folder                        # --prefix e-1-paroptqn
                )
                case_cmds.append(cmd)
            else:
                n_exist_case += 1

    with open(args.output, mode='w') as f:
        for row in case_cmds:
            f.write(row + '\n')

    # Print out summary
    n_planned_cases = len(args.optimizer)*len(cases)
    n_new_cases = len(case_cmds)
    print("------ Summary ------")
    print("Total number of cases needed:{:4d}".format(n_planned_cases))
    print("Number of existing cases:    {:4d}".format(n_exist_case))
    print("Number of new cases created: {:4d}".format(n_new_cases))

    if(n_new_cases + n_exist_case != n_planned_cases):
        raise RuntimeError('Case numbers don\'t match!')

    # Generate pbs
    with open('eig_cases.pbs', 'w') as f:
        f.write('#PBS -N eig_n{:d}\n'.format(n_new_cases))
        f.write('#PBS -A GT-gkennedy9-CODA20\n')
        f.write('#PBS -l nodes=1:ppn=24\n')
        f.write('#PBS -l pmem=6gb\n')
        f.write('#PBS -l walltime={:d}:00:00\n'.format(args.walltime))
        f.write('#PBS -j oe\n')
        f.write('#PBS -o cases-n{:d}.out\n'.format(n_new_cases))
        f.write('#PBS -m a\n')
        f.write('#PBS -M yfu97@gatech.edu\n')
        f.write('#PBS -t 1-{:d}\n'.format(n_new_cases))
        f.write('\n')
        f.write('date\n')
        f.write('\n')
        f.write('cd $PBS_O_WORKDIR\n')
        f.write('\n')
        f.write('#Get ith line of target file and run\n')
        f.write('export casecmd=`cat eig_cases.txt | head -$PBS_ARRAYID | tail -1`\n')
        f.write('\n')
        f.write('mpirun -np 24 $casecmd\n')
        f.write('\n')
        f.write('date\n')
        f.write('\n')

