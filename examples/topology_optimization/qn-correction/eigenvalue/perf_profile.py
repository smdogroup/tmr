import numpy as np
import argparse
import os
import re
import pickle
import matplotlib.pyplot as plt
from pprint import pprint
from csv import DictWriter

PERF_INF = 1e20
optimizers = ['paropt', 'paroptqn', 'snopt', 'ipopt']
colors = {
    'paropt': '#00876C',
    'paroptqn': '#BC5090',
    'ipopt': '#7A4EFE',
    'snopt': '#2e2e2e',
}
legends = {
    'paropt': r'ParOpt filterSQP',
    'paroptqn': r'ParOpt filterSQP w/ qn',
    'ipopt': r'IPOPT',
    'snopt': r'SNOPT',
}


def getDirs(result_folder):
    """
    get all folder names with the following pattern:

        e-number-name

    for example:

        e-12-paroptqn

    Args:
        result_folder (str): the folder conatining all `e-number-optimizer' subfolders

    Return:
        dirs (list): list of dir names
    """

    dirs = os.listdir(result_folder)
    r = re.compile(r"e-\d+-.+")
    dirs = list(filter(r.match, dirs))

    return dirs

def createDicStruct(dirs):
    """
    Create the main dictionary to get data stored
    in result pkl files

    Args:
        dirs (list): list of dir names

    Return:
        physics (dict): an empty dictionary structure,  note that its
                        length is the number of physical problems
    """

    physics = dict()
    for d in dirs:
        _, num, omz = d.split('-')
        if num not in physics.keys():
            physics[num] = dict()
        physics[num][omz] = dict()
        physics[num][omz]['obj'] = None
        physics[num][omz]['obj_normed'] = None
        physics[num][omz]['infeas'] = None
        physics[num][omz]['discreteness'] = None
        physics[num]['best_obj'] = PERF_INF

    return physics

def populateDic(dirs, n_mesh_refine, result_folder, physics):
    """
    Populate the dictioanry structure

    Args:
        dirs (list): list of dir names
        n_mesh_refine (int): number of mesh refinements, this decides which
                             pickle contains the final result
        result_folder (str): the folder conatining all `e-number-optimizer' subfolders
        physics (dict): the dictionary structure to populate

    Return:
        n_pkls (int): number of existing pickles
    """

    # We keep tracking number of existing pickles
    n_pkls = 0

    for d in dirs:
        _, num, omz = d.split('-')
        pklname = 'output_refine{:d}.pkl'.format(n_mesh_refine-1)
        pklpath = os.path.join(args.result_folder, d, pklname)

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

def normalizeDict(physics, infeas_tol):
    """
    normalize the objective against the best

    Args:
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
                        physics[num][omz]['obj_normed'] = -PERF_INF
                        print('[Info] Infeasibile case detected, No:{:>5s}, ' \
                              'optimizer:{:>10s}, infeas:{:>20.10e}'.format(
                               num, omz, physics[num][omz]['infeas']))
    return n_feas

def genProfileData(optimizers, physics):
    """
    Prepare raw data for the profile

    Args:
        optimizers (list): list of optimizers used
        physics (dict): the dictionary structure

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
        profile_data[omz]['obj_normed'].sort(reverse=True)
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

def plotObjProfile(profile_data, obj_bound, optimizers):
    """
    Plot objective profile

    Args:
        profile_data (dict): profile data
        obj_bound (float): smallest objective to be included in the profile
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
        if x[-1] > obj_bound:
            x.append(obj_bound)
            y.append(y[-1])
        ax.step(x, y, label=legends[omz], color=colors[omz])

    ax.set_xlim(1.0 + 0.05*(1.0 - obj_bound), obj_bound)
    ax.set_xlabel('Normalized objective (eigenvalue)')
    ax.set_ylabel('Fraction of cases')
    fig.legend(loc='center right')

    return fig, ax

def plotDiscreteProfile(profile_data, optimizers):
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

def saveTable(physics, result_folder, csv='physics.csv'):
    """
    Save physics to a csv in result folder

    Args:
        physics (dict): main dictionary
        csv (str): csv name
    """

    with open(os.path.join(result_folder, csv), 'w') as f:
        fieldnames = ['no', 'paropt', 'paroptqn', 'snopt', 'ipopt']
        writer = DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for num in physics:
            writer.writerow({
                'no': num,
                'paropt': physics[num]['paropt']['obj'],
                'paroptqn': physics[num]['paroptqn']['obj'],
                'snopt': physics[num]['snopt']['obj'],
                'ipopt': physics[num]['ipopt']['obj'],
            })

    return

if __name__ == '__main__':

    # Argument parser
    p = argparse.ArgumentParser()
    p.add_argument('--result-folder', type=str, default='.')
    p.add_argument('--infeas-tol', type=float, default=1e-6)
    p.add_argument('--n-mesh-refine', type=int, default=1)
    p.add_argument('--obj-bound', type=float, default=0.5)
    p.add_argument('--obj-profile-name', type=str, default='profile-obj.pdf')
    p.add_argument('--dis-profile-name', type=str, default='profile-discrete.pdf')
    args = p.parse_args()

    # Get list of case dirs by matching folder name
    dirs = getDirs(args.result_folder)

    # Create dictionary structure
    physics = createDicStruct(dirs)

    # Populate dictionary
    n_pkls = populateDic(dirs, args.n_mesh_refine, args.result_folder, physics)

    # Normalize objective
    n_feas = normalizeDict(physics, args.infeas_tol)

    # Generate profile data
    profile_data, n_best = genProfileData(optimizers, physics)

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
    fig, ax = plotObjProfile(profile_data, args.obj_bound, optimizers)
    fig.savefig(os.path.join(args.result_folder, args.obj_profile_name))

    # Plot discreteness profile
    fig, ax = plotDiscreteProfile(profile_data, optimizers)
    fig.savefig(os.path.join(args.result_folder, args.dis_profile_name))

    # Save physics as csv
    saveTable(physics, args.result_folder)
