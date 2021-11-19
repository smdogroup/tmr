import numpy as np
import argparse
import os
import re
import pickle
import matplotlib.pyplot as plt
from pprint import pprint
from csv import DictWriter

PERF_INF = 1e20
colors = {
    'paropt': '#00876C',
    'paroptqn': '#BC5090',
    'ipopt': '#7A4EFE',
    'snopt': '#2e2e2e',
    'mma': '#FFA600'
}

def getDirs(result_folders):
    """
    get all folder names with the following pattern:

        number-name

    for example:

        12-paroptqn

    Args:
        result_folders (list): the list of folders conatining all `number-optimizer' subfolders

    Return:
        dirslist (list): list of dir names lists
    """

    dirslist = []
    for f in result_folders:
        dirs = os.listdir(f)
        r = re.compile(r"\d+-.+")
        dirs = list(filter(r.match, dirs))
        sort_key = lambda text : int(text.split('-')[0])
        dirs.sort(key=sort_key)
        dirslist.append(dirs)

    return dirslist

def createDicStruct(dirslist):
    """
    Create the main dictionary to get data stored
    in result pkl files

    Args:
        dirslist (list): list of dir names

    Return:
        physics (dict): an empty dictionary structure,  note that its
                        length is the number of physical problems
    """

    physics = dict()
    for dirs in dirslist:
        for d in dirs:
            num, omz = d.split('-')
            if num not in physics.keys():
                physics[num] = dict()
            physics[num][omz] = dict()
            physics[num][omz]['obj'] = None
            physics[num][omz]['obj_normed'] = None
            physics[num][omz]['infeas'] = None
            physics[num][omz]['discreteness'] = None
            physics[num]['best_obj'] = PERF_INF

    return physics

def populateDic(dirslist, n_mesh_refine, result_folders, physics):
    """
    Populate the dictioanry structure

    Args:
        dirslist (list): list of dir names lists
        n_mesh_refine (int): number of mesh refinements, this decides which
                             pickle contains the final result
        result_folders (list): the list of folders conatining all `number-optimizer' subfolders
        physics (dict): the dictionary structure to populate

    Return:
        n_pkls (int): number of existing pickles
    """

    # We keep tracking number of existing pickles
    n_pkls = 0

    for i, dirs in enumerate(dirslist):
        for d in dirs:
            num, omz = d.split('-')
            pklname = 'output_refine{:d}.pkl'.format(n_mesh_refine-1)
            pklpath = os.path.join(result_folders[i], d, pklname)

            # Success case
            try:
                with open(pklpath, 'rb') as f:
                    pkldict = pickle.load(f)
                    physics[num][omz]['obj'] = pkldict['obj']
                    physics[num][omz]['infeas'] = pkldict['infeas']
                    physics[num][omz]['discreteness'] = pkldict['discreteness']

                    if pkldict['obj'] < physics[num]['best_obj']:
                        physics[num]['best_obj'] = pkldict['obj']

                    n_pkls += 1

            # Fail case, no pkl generated
            except:
                print('[Info] {:} doesn\'t exist!'.format(pklpath))
                physics[num][omz]['obj'] = PERF_INF
                physics[num][omz]['infeas'] = PERF_INF
                physics[num][omz]['discreteness'] = PERF_INF

    return n_pkls

def normalizeDict(problem, physics, infeas_tol):
    """
    normalize the objective against the best

    Args:
        problem (str): 'eig' or 'comp'
        physics (dict): the dictionary structure
        infeas_tol (float): maximum acceptable infeasibility

    Return:
        n_feas (int): number of feasible cases
    """

    # Keep tracking the feasible cases
    n_feas = 0

    # Normalize obj against best
    for num in physics.keys():
        for omz in physics[num]:
            if omz != 'best_obj':
                if physics[num][omz]['obj'] is not None:
                    if physics[num][omz]['infeas'] < infeas_tol:
                        physics[num][omz]['obj_normed'] = \
                            physics[num][omz]['obj'] / physics[num]['best_obj']
                        n_feas += 1
                    else:
                        if problem == 'eig':
                            physics[num][omz]['obj_normed'] = -PERF_INF
                        else:
                            physics[num][omz]['obj_normed'] = PERF_INF
                        physics[num][omz]['discreteness'] = PERF_INF
                        if physics[num][omz]['infeas'] < PERF_INF:
                            print('[Info] Infeasibile case detected, No:{:>5s}, ' \
                                'optimizer:{:>10s}, infeas:{:>20.10e}'.format(
                                num, omz, physics[num][omz]['infeas']))
    return n_feas

def genProfileData(optimizers, physics, problem):
    """
    Prepare raw data for the profile

    Args:
        optimizers (list): list of optimizers used
        physics (dict): the dictionary structure
        problem (str): 'eig' or 'comp'

    Return:
        profile_data (dict): profile data dictionary
        n_best (int): number of winners, this should be equal to number of physical problems
    """

    # Dimensinless objective for each optimizer
    profile_data = {}
    n_physics = len(physics)
    for omz in optimizers:
        profile_data[omz] = dict()

        # Prepare objective data
        profile_data[omz]['obj_normed'] = [physics[num][omz]['obj_normed'] for num in physics]
        if problem == 'eig':
            profile_data[omz]['obj_normed'].sort(reverse=True)
        else:
            profile_data[omz]['obj_normed'].sort(reverse=False)
        profile_data[omz]['obj_percentile'] = [(index+1)/n_physics for index, _ in enumerate(profile_data[omz]['obj_normed'])]

        # Prepare discreteness data
        profile_data[omz]['dis'] = [physics[num][omz]['discreteness'] for num in physics]
        profile_data[omz]['dis'].sort(reverse=False)
        profile_data[omz]['dis_percentile'] = [(index+1)/n_physics for index, _ in enumerate(profile_data[omz]['dis'])]


    # Check if number of 1.0 equals number of physical problem
    n_best = 0
    for omz in optimizers:
        for obj in profile_data[omz]['obj_normed']:
            if obj == 1.0:
                n_best += 1

    return profile_data, n_best

def plotObjProfile(problem, profile_data, eig_bound, comp_bound, optimizers, legends):
    """
    Plot objective profile

    Args:
        problem (str): 'eig' or 'comp'
        profile_data (dict): profile data
        eig_bound (float): smallest objective to be included in the profile
        comp_bound (float): largest objective to be included in the profile
        optimizers (list): list of optimizers used

    Return:
        fig (matplotlib.Figure)
        ax (matplotlib.axes.Axes)
    """

    fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8))
    for omz in optimizers:
        # Append an artificial entry to the list so we have
        # nice-looking profile at the end
        x = profile_data[omz]['obj_normed'].copy()
        y = profile_data[omz]['obj_percentile'].copy()
        if problem == 'eig':
            if x[-1] > eig_bound:
                x.append(eig_bound)
                y.append(y[-1])
        else:
            if x[-1] < comp_bound:
                x.append(comp_bound)
                y.append(y[-1])

        ax.step(x, y, label=legends[omz], color=colors[omz])

    if problem == 'eig':
        ax.set_xlim(1.0 + 0.05*(1.0 - eig_bound), eig_bound)
    else:
        ax.set_xlim(1.0 - 0.05*(comp_bound - 1.0), comp_bound)

    if problem == 'eig':
        ax.set_xlabel('Normalized objective (eigenvalue)')
    elif problem == 'comp':
        ax.set_xlabel('Normalized objective (compliance)')
    elif problem == 'freq':
        ax.set_xlabel('Normalized objective (mass)')
    ax.set_ylabel('Fraction of cases')
    fig.legend(loc='center right')

    return fig, ax

def plotDiscreteProfile(profile_data, optimizers, legends):
    """
    Plot discreteness profile

    Args:
        profile_data (dict): profile data
        optimizers (list): list of optimizers used

    Return:
        fig (matplotlib.Figure)
        ax (matplotlib.axes.Axes)
    """

    fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8))
    for omz in optimizers:
        # Append an artificial entry to the list so we have
        # nice-looking profile at the end
        x = profile_data[omz]['dis'].copy()
        y = profile_data[omz]['dis_percentile'].copy()
        x.append(0.25)
        y.append(y[-1])
        ax.step(x, y, label=legends[omz], color=colors[omz])

    ax.set_xlim(-0.01, 0.25)
    ax.set_xlabel('Averaged discreteness')
    ax.set_ylabel('Fraction of cases')
    fig.legend(loc='center right')

    return fig, ax

def omzFields(omz):
    return [omz+'-obj', omz+'-infeas', omz+'-discrete']

def saveTable(physics, result_folder, optimizers, csv='physics.csv'):
    """
    Save physics to a csv in result folder

    Args:
        physics (dict): main dictionary
        csv (str): csv name
    """

    fieldnames = ['no']
    for omz in optimizers:
        fieldnames.extend(omzFields(omz))

    with open(os.path.join(result_folder, csv), 'w') as f:
        writer = DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for num in physics:
            write_row_dict = {'no': num}
            for omz in optimizers:
                write_row_dict[omz+'-obj'] = physics[num][omz]['obj']
                write_row_dict[omz+'-infeas'] = physics[num][omz]['infeas']
                write_row_dict[omz+'-discrete'] = physics[num][omz]['discreteness']

            writer.writerow(write_row_dict)

    return

if __name__ == '__main__':

    # Argument parser
    p = argparse.ArgumentParser()
    p.add_argument('--result-folder', type=str, nargs='*', default=['.'])
    p.add_argument('--problem', type=str, default='eig', choices=['eig', 'comp', 'freq'])
    p.add_argument('--infeas-tol', type=float, default=1e-6)
    p.add_argument('--optimizers', type=str, nargs='*',
        default=['paropt', 'paroptqn', 'snopt', 'ipopt', 'mma'],
        choices=['paropt', 'paroptqn', 'snopt', 'ipopt', 'mma'])
    p.add_argument('--n-mesh-refine', type=int, default=1)
    p.add_argument('--eig-bound', type=float, default=0.5)
    p.add_argument('--comp-bound', type=float, default=1.5)
    p.add_argument('--paropt-type', type=str, default='sl1QP',
        choices=['sl1QP', 'filterSQP'])
    args = p.parse_args()

    legends = {
        'paropt': r'ParOpt {:s}'.format(args.paropt_type),
        'paroptqn': r'ParOpt {:s} w/ qn'.format(args.paropt_type),
        'ipopt': r'IPOPT',
        'snopt': r'SNOPT',
        'mma': r'MMA'
    }
    legends = {key: legends[key] for key in args.optimizers}

    # Get list of case dirs by matching folder name
    dirslist = getDirs(args.result_folder)

    # Create dictionary structure
    physics = createDicStruct(dirslist)

    # Populate dictionary
    n_pkls = populateDic(dirslist, args.n_mesh_refine, args.result_folder, physics)

    # Normalize objective
    n_feas = normalizeDict(args.problem, physics, args.infeas_tol)

    # Generate profile data
    optimizers = args.optimizers
    profile_data, n_best = genProfileData(optimizers, physics, args.problem)

    # Print summary
    n_physics = len(physics)
    n_planned = len(physics)*len(optimizers)
    print('\n')
    print('-----------SUMMARY-----------')
    print('number of physics:         {:d}'.format(n_physics))
    print('number of best cases       {:d} (should == {:d})'.format(n_best, n_physics))
    print('number of planned cases:   {:d}'.format(n_planned))
    print('number of existing pkls:   {:d} ({:.2f}%)'.format(n_pkls, 100*n_pkls/n_planned))
    print('number of feasible cases:  {:d} ({:.2f}%)'.format(n_feas, 100*n_feas/n_planned))

    # Set up plotting environment
    mpl_style_path = os.path.dirname(os.path.realpath(__file__)) + '/paper.mplstyle'
    plt.style.use(mpl_style_path)

    # Plot objective profile
    fig, ax = plotObjProfile(args.problem, profile_data, args.eig_bound,
        args.comp_bound, optimizers, legends)
    if len(args.result_folder) == 1:
        fig.savefig(os.path.join(args.result_folder[0], 
            "profile_obj_infeas_{:.0e}.pdf".format(args.infeas_tol)), dpi=800)
    else:
        fig.savefig("profile_obj_infeas_{:.0e}.pdf".format(args.infeas_tol), dpi=800)

    # Plot discreteness profile
    fig, ax = plotDiscreteProfile(profile_data, optimizers, legends)
    if len(args.result_folder) == 1:
        fig.savefig(os.path.join(args.result_folder[0],
            "profile_dis_infeas_{:.0e}.pdf".format(args.infeas_tol)), dpi=800)
    else:
        fig.savefig("profile_dis_infeas_{:.0e}.pdf".format(args.infeas_tol), dpi=800)

    # Save physics as csv
    if len(args.result_folder) == 1:
        saveTable(physics, args.result_folder[0], args.optimizers)
    else:
        saveTable(physics, '.', args.optimizers)
