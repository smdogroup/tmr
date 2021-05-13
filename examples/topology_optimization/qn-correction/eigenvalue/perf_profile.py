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

def plotObj(data, obj_bound):
    fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8))
    for omz in optimizers:
        ax.step(data[omz]['obj_normed'], data[omz]['percentile'],
                label=legends[omz], color=colors[omz])

    ax.set_xlim(1.0 + 0.05*(1.0 - obj_bound), obj_bound)
    ax.set_xlabel('Normalized objective (eigenvalue)')
    ax.set_ylabel('Fraction of cases')
    fig.legend(loc='center right')
    fig.savefig(os.path.join(args.result_folder, 'obj.pdf'))
    return

def saveTable(physics, csv='eig_data.csv'):
    with open(os.path.join(args.result_folder,csv), 'w') as f:
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
    p.add_argument('--result-folder', type=str, default='./pkls')
    p.add_argument('--infeas-tol', type=float, default=1e-6)
    p.add_argument('--n-mesh-refine', type=int, default=1)
    p.add_argument('--obj-bound', type=float, default=0.5)
    p.add_argument('--dis-bound', type=float, default=0.25)
    args = p.parse_args()

    # Get list of case dirs by matching folder name
    dirs = os.listdir(args.result_folder)
    r = re.compile(r"e-\d+-.+")
    dirs = list(filter(r.match, dirs))

    # Create dictionary structure
    physics = dict()
    for d in dirs:
        e, num, omz = d.split('-')
        if num not in physics.keys():
            physics[num] = dict()
        physics[num][omz] = dict()
        physics[num][omz]['obj'] = None
        physics[num][omz]['obj_normed'] = None
        physics[num][omz]['infeas'] = None
        physics[num][omz]['discreteness'] = None
        physics[num]['best_obj'] = PERF_INF

    # Keep track numbers
    n_physics = len(physics)
    n_planned = len(physics)*len(optimizers)
    n_pkls = 0
    n_feas = 0

    # Populate dictionary
    for d in dirs:
        e, num, omz = d.split('-')
        pklname = 'output_refine{:d}.pkl'.format(args.n_mesh_refine-1)
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

    # Normalize obj against best
    for num in physics.keys():
        for omz in physics[num]:
            if omz != 'best_obj':
                if physics[num][omz]['obj'] is not None:
                    if physics[num][omz]['infeas'] < args.infeas_tol:
                        physics[num][omz]['obj_normed'] = \
                            physics[num][omz]['obj'] / physics[num]['best_obj']
                        n_feas += 1
                    else:
                        physics[num][omz]['obj_normed'] = -PERF_INF
                        print('[Info] Infeasibile case detected, No:{:>5s}, ' \
                              'optimizer:{:>10s}, infeas:{:>20.10e}'.format(
                               num, omz, physics[num][omz]['infeas']))

    # Dimensinless objective for each optimizer
    data = {}
    for omz in optimizers:
        data[omz] = dict()
        data[omz]['obj_normed'] = [physics[num][omz]['obj_normed'] for num in physics]
        data[omz]['obj_normed'].sort(reverse=True)
        data[omz]['percentile'] = [(index+1)/n_physics for index, _ in enumerate(data[omz]['obj_normed'])]


    # Check if number of 1.0 equals number of physical problem
    n_best = 0
    for omz in optimizers:
        for obj in data[omz]['obj_normed']:
            if obj == 1.0:
                n_best += 1

    # Print summary
    print('\n')
    print('-----------SUMMARY-----------')
    print('number of physics:         {:d}'.format(n_physics))
    print('number of best cases       {:d}'.format(n_best))
    print('number of planned cases:   {:d}'.format(n_planned))
    print('number of existing pkls:   {:d} ({:.2f}%)'.format(n_pkls, 100*n_pkls/n_planned))
    print('number of feasible cases:  {:d} ({:.2f}%)'.format(n_feas, 100*n_feas/n_planned))

    # Set up plotting environment
    mpl_style_path = os.path.dirname(os.path.realpath(__file__)) + '/paper.mplstyle'
    plt.style.use(mpl_style_path)

    plotObj(data, args.obj_bound)
    saveTable(physics)









