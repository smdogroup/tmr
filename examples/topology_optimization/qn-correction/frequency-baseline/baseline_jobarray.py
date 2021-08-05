from csv import DictWriter
from typing import ClassVar
from baseline import run_baseline_case
import argparse

import sys
sys.path.append('../eigenvalue')
from job_array import readCSV, dic2str

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

def writeToCSV(csvfile, ls, fieldnames):
    # Create csv dictionary reader
    with open(csvfile, mode='w', encoding='utf-8-sig') as f:
        writer = DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in ls:
            writer.writerow(row)
    return

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    p.add_argument("--csv", type=str, default='baselines.csv')
    p.add_argument("--write_result_to_csv", action='store_true')
    args = p.parse_args()

    physics = readCSV(args.csv, args.start, args.end)

    results = []
    for case_dict in physics:
        for k in case_dict:
            if isint(case_dict[k]):
                case_dict[k] = int(case_dict[k])
            elif isfloat(case_dict[k]):
                case_dict[k] = float(case_dict[k])
        res_dict = case_dict.copy()
        case_dict.pop('no')
        res = run_baseline_case(**case_dict)
        res_dict['min eig'] = res['min eig']
        results.append(res_dict)

    fieldnames = [s for s in results[0]]
    writeToCSV('baseline_results.csv', results, fieldnames)






