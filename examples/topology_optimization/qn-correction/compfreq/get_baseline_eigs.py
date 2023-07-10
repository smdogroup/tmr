"""
This script extracts baseline frequencies from pkl files and
store the result in csv
"""

import re
import os
import pickle
from csv import DictWriter
import argparse


class Extractor:
    def __init__(self, result_folder, start, end):
        self.result_folder = result_folder
        self.start = start
        self.end = end
        return

    def getBLEigCSV(self, csv="baseline_frq.csv"):
        cols = ["no", "paropt", "paroptqn", "snopt", "ipopt"]
        optimizers = ["paropt", "paroptqn", "snopt", "ipopt"]
        with open(os.path.join(self.result_folder, csv), "w") as csvf:
            writer = DictWriter(csvf, fieldnames=cols)
            writer.writeheader()
            for i in range(self.start, self.end + 1):
                row_dict = {"no": i}
                for omz in optimizers:
                    try:
                        with open(
                            os.path.join(
                                self.result_folder,
                                "e-{:d}-{:s}".format(i, omz),
                                "output_refine0.pkl",
                            ),
                            "rb",
                        ) as pklf:
                            eig = pickle.load(pklf)["min_eig"]
                    except:
                        eig = "FAIL_LOAD_PKL"
                    row_dict[omz] = eig
                writer.writerow(row_dict)
        return


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("result_folder", type=str)
    p.add_argument("start", type=int)
    p.add_argument("end", type=int)
    args = p.parse_args()

    extractor = Extractor(args.result_folder, args.start, args.end)
    extractor.getBLEigCSV()
