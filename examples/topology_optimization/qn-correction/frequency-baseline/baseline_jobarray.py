from csv import DictWriter, DictReader
from typing import ClassVar
from baseline import run_baseline_case
import argparse
from mpi4py import MPI

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

def isfloat(x):
    try:
        a = float(x)
    except (TypeError, ValueError):
        return False
    else:
        return True

def isint(x):
    try:
        a = float(x)
        b = int(a)
    except (TypeError, ValueError):
        return False
    else:
        return a == b

def writeLineToCSV(csvfile, line, fieldnames, write_header=False):
    # Create csv dictionary reader
    with open(csvfile, mode='a', encoding='utf-8-sig') as f:
        writer = DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(line)
    return

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    p.add_argument("--csv", type=str, default='cases.csv')
    p.add_argument("--write_result_to_csv", action='store_true')
    args = p.parse_args()

    comm = MPI.COMM_WORLD
    physics = readCSV(args.csv, args.start, args.end)

    index = 0
    for case_dict in physics:
        for k in case_dict:
            if isint(case_dict[k]):
                case_dict[k] = int(case_dict[k])
            elif isfloat(case_dict[k]):
                case_dict[k] = float(case_dict[k])
        res_dict = case_dict.copy()
        case_dict.pop('no')
        try:
            res = run_baseline_case(**case_dict)
            res_dict['min eig'] = res['min eig']
            res_dict['sol time'] = '{:.2f} s'.format(res['sol time'])
        except Exception as e:
            res_dict['min eig'] = 'case failed!'

        fieldnames = [s for s in res_dict]

        write_header = False
        if index == 0:
            write_header = True

        if comm.rank == 0:
            writeLineToCSV('results.csv', res_dict, fieldnames, write_header)

        index += 1

