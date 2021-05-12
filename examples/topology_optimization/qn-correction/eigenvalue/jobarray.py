"""
Generate a bash script to run multiple serial cases concurrently using gnu-parallel
"""

import os
import sys
import argparse
import numpy as np

def optimizer(omz):
    if omz == 'paropt':
        return 'paropt-pyoptsparse'
    elif omz == 'paroptqn':
        return 'paropt-pyoptsparse --qn-correction'
    else:
        return omz

def domname(dom):
    if dom == 'lbracket':
        return 'lbr'
    elif dom == 'cantilever':
        return 'can'
    elif dom == 'mbb':
        return 'mbb'
    elif dom == 'michell':
        return 'mic'

def list2str(ls):
    '''
    Convert a list to a single string divided by space
    example:
    ls = ['aaa', 'bbb', 2]
    list2str(ls) = 'aaa bbb 2'
    '''
    return ' '.join([str(e) for e in ls])

if __name__ == '__main__':

    # Create the argument parser
    p = argparse.ArgumentParser()

    # Analysis list parameters
    p.add_argument('--domain', type=str, nargs='*',
        default=['cantilever', 'michell', 'mbb', 'lbracket'],
        choices=['cantilever', 'michell', 'mbb', 'lbracket'])
    p.add_argument('--AR', type=float, nargs='*', default=[1.0, 2.0, 3.0, 4.0, 5.0])
    p.add_argument('--ratio', type=float, nargs='*', default=[0.3, 0.4, 0.5, 0.6, 0.7])

    # Optimization list parameters
    p.add_argument('--vol-frac', type=float, nargs='*', default=[0.4])
    p.add_argument('--optimizer', type=str, nargs='*',
        default=['paropt', 'paroptqn', 'snopt', 'ipopt'],
        choices=['paropt', 'paroptqn', 'snopt', 'ipopt'])

    # Analysis scalar parameters
    # p.add_argument('--r0-frac', type=float, default=0.05)
    # p.add_argument('--mg-levels', type=int, default=4)
    p.add_argument('--qval', type=float, default=5.0)
    p.add_argument('--max-jd-size', type=int, default=100)
    p.add_argument('--max-gmres-size', type=int, default=30)

    # Optimization scalar parameters
    p.add_argument('--n-mesh-refine', type=int, default=1)
    p.add_argument('--max-iter', type=int, default=100)
    p.add_argument('--qn-subspace', type=int, default=2)

    # Executable directory
    p.add_argument('--exec', type=str, default='eig-max.py')
    p.add_argument('--create-models-exec', type=str, default='create_models.py')

    # PBS settings
    # p.add_argument('--partition', type=int, default=1)
    p.add_argument('--jobname', type=str, default='serial')
    p.add_argument('--walltime', type=int, default=24)
    p.add_argument('--pmem', type=int, default=6)

    # Parse arguments
    args = p.parse_args()

    # Save the command and arguments that executed this script
    cmd = 'python ' + ' '.join(sys.argv)

    # Create egads models
    if 'lbracket' not in args.domain:
        egads_domains = ['cantilever']
    else:
        egads_domains = ['cantilever', 'lbracket']
    create_models_exec = 'python {:s} --domain {:s} --AR {:s} --ratio {:s}'.format(
        args.create_models_exec, list2str(egads_domains), list2str(args.AR), list2str(args.ratio))
    print('[generate_scripts.py] Creating egads models...')
    os.system(create_models_exec)


    # Count total number of cases, notice that --AR option is only valid
    # for cantilever, michell and mbb, and --ratio option is only valid for lbracket
    n_domain = len(args.domain)
    n_AR     = len(args.AR)
    n_ratio  = len(args.ratio)
    n_opt    = len(args.optimizer)
    n_vol    = len(args.vol_frac)
    if 'lbracket' not in args.domain:
        n_cases = n_domain * n_AR * n_opt * n_vol
    else:
        n_cases = (n_domain-1)*n_AR + 1*n_ratio
        n_cases *= n_opt*n_vol

    print('[generate_scripts.py] {:d} cases needed...'.format(n_cases))

    # Create a string contain all scalar parameters
    commons = "--r0-frac {:f} --mg-levels {:d} --qval {:f} --max-jd-size {:d} " \
              "--max-gmres-size {:d} --n-mesh-refine {:d} --max-iter {:d} " \
              "--qn-subspace {:d}".format(args.r0_frac, args.mg_levels, args.qval,
              args.max_jd_size, args.max_gmres_size, args.n_mesh_refine, args.max_iter,
              args.qn_subspace)

    # Create a list containting all case entries
    cases = []
    case_num = 0
    case_exist = 0
    for dom in args.domain:
        for vol in args.vol_frac:
            for omz in args.optimizer:

                # If this is an lbracket domain
                if dom == 'lbracket':
                    for rat in args.ratio:
                        prefix = 'lbr-r{:.1f}-v{:.1f}-{:s}'.format(rat, vol, omz)
                        if not os.path.isdir(prefix):
                            os.mkdir(prefix)
                        if os.path.isfile(os.path.join(prefix, 'output_refine0.pkl')):
                            case_exist += 1
                        else:
                            options = "--domain {:s} --AR 1.0 --ratio {:f} --vol-frac {:f} " \
                                    "--optimizer {:s} --htarget {:f} --prefix {:s}".format(
                                    dom, rat, vol, optimizer(omz), 0.5*rat, prefix)
                            case_num += 1
                            cases.append(" ".join(["echo \"[{:4d}/{:4d}] running case: {:s}\" && python".format(
                                        case_num, n_cases, prefix), args.exec, options, commons, ">",
                                        os.path.join(prefix, "{:s}.out".format(prefix))]))

                # Otherwise, this is a cantilever-type domain
                else:
                    for ar in args.AR:
                        prefix = '{:s}-a{:.1f}-v{:.1f}-{:s}'.format(domname(dom), ar, vol, omz)
                        if not os.path.isdir(prefix):
                            os.mkdir(prefix)
                        if os.path.isfile(os.path.join(prefix, 'output_refine0.pkl')):
                            case_exist += 1
                        else:
                            options = "--domain {:s} --AR {:f} --ratio 1.0 --vol-frac {:f} " \
                                    "--optimizer {:s} --htarget 0.5 --prefix {:s}".format(
                                    dom, ar, vol, optimizer(omz), prefix)
                            case_num += 1
                            cases.append(" ".join(["echo \"[{:4d}/{:4d}] running case: {:s}\" && python".format(
                                        case_num, n_cases, prefix), args.exec, options, commons, ">",
                                        os.path.join(prefix, "{:s}.out".format(prefix))]))

    # Write lines to shell script file
    sh_name = 'cases-n{:d}.sh'.format(case_num)
    if os.path.isfile(sh_name):
        print ("[generate_scripts.py] Warning: File {:s} exists! Replacing it...".format(sh_name))

    with open(sh_name, 'w') as f:
        f.write("# Command that generates this file:\n# {:s}\n".format(cmd))
        for entry in cases:
            f.write(entry + '\n')

    # Print out summary
    print("------ Summary ------")
    print("Total number of cases needed:{:4d}".format(n_cases))
    print("Number of existing cases:    {:4d}".format(case_exist))
    print("Number of new cases created: {:4d}".format(case_num))

    exit(0)

    # Compute number of cases in each partition
    part_nums = np.zeros(args.partition)
    part_nums[:] = case_num // args.partition
    for i in range(case_num % args.partition):
        part_nums[i] += 1
    print("Total number in partitions:  {:4d}".format(int(np.sum(part_nums))))

    # Partition the script and create PBS scripts
    global_index = 0
    for pi in range(args.partition):

        # Create partition shell scripts
        part_name = 'cases-n{:d}-part{:d}.sh'.format(case_num, pi+1)
        with open(part_name, 'w') as f:
            for local_index in range(int(part_nums[pi])):
                f.write(cases[global_index] + '\n')
                global_index += 1

        # Create corresponding pbs scripts
        pbs_name = 'cases-n{:d}-part{:d}.pbs'.format(case_num, pi+1)
        with open(pbs_name, 'w') as f:
            f.write('#PBS -N {:s}-{:d}\n'.format(args.jobname, pi+1))
            f.write('#PBS -A GT-gkennedy9-CODA20\n')
            f.write('#PBS -l nodes=1:ppn=24\n')
            f.write('#PBS -l pmem={:d}gb\n'.format(args.pmem))
            f.write('#PBS -l walltime={:d}:00:00\n'.format(args.walltime))
            f.write('#PBS -j oe\n')
            f.write('#PBS -o cases-n{:d}-part{:d}.out\n'.format(case_num, pi+1))
            f.write('#PBS -m abe\n')
            f.write('\n')
            f.write('cd $PBS_O_WORKDIR\n')
            f.write('parallel -j 24 :::: cases-n{:d}-part{:d}.sh\n'.format(case_num, pi+1))

    # Generate submit
    with open('submit', 'w') as f:
        for pi in range(args.partition):
            f.write('qsub cases-n{:d}-part{:d}.pbs\n'.format(case_num, pi+1))

    os.system('chmod +x submit')












