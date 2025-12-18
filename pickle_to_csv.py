import os, pickle
import pandas as pd


def main(dirname="./small_models"):
    filenames = os.listdir(dirname)
    for filename in filenames:
        filepath = os.path.join(dirname, filename)
        if os.path.isfile(filepath) and filename.endswith(".pickle"):
            print(filepath)
            df = pickle.load(open(filepath, "rb"))
            df = pd.DataFrame(df).T
            df.to_csv(filename.replace(".pickle", ".csv"), index=True)


if "__main__" == __name__:
    import fire

    fire.Fire(main)
