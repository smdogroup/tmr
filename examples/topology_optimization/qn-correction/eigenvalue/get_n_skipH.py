import pickle
import argparse
import re
import os


def getDirs(result_folder):
    """
    get all folder names with the following pattern:

        e-number-name

    for example:

        e-12-paroptqn

    Args:
        result_folder (list): folder conatining all `e-number-optimizer' subfolders

    Return:
        dirs (list): list of dir names
    """

    dirs = os.listdir(result_folder)
    r = re.compile(r"e-\d+-.+")
    dirs = list(filter(r.match, dirs))
    dirs.sort()

    return dirs


if __name__ == "__main__":
    # Argument parser
    p = argparse.ArgumentParser()
    p.add_argument("--result-folder", type=str, default=["."])
    args = p.parse_args()

    dirs = getDirs(args.result_folder)

    for d in dirs:
        if "paroptqn" in d:
            pklname = os.path.join(args.result_folder, d, "output_refine0.pkl")

            with open(pklname, "rb") as f:
                pkldict = pickle.load(f)
                nskipH = pkldict["n_skipH"]
                print("{:<20s}{:<20s}{:<20d}".format(d, "nskipH", nskipH))
