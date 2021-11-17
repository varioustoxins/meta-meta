
import pathlib
import os
import csv
import sys

class Match:
    def __init__(self):
        self.shifts = {}
        self.names = {}

    def load_shifts(self,path):

        path =pathlib.Path(path,'id_shifts.csv')
        with open(path, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='"')
            for row in reader:
                self.shifts[(row[0],row[1])] = [float(elem) for elem in row[2:]]

    def load_names(self,path):

        path =pathlib.Path(path,'id_name.csv')
        with open(path, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='"')
            for row in reader:
                self.names[row[0]] = row[1]

    def load(self):
        cwd =  os.getcwd()
        for directory in ['hmdb_nmr_spectra','mmcd_nmr_spectra','bmrb_nmr_spectra']:
            path = (pathlib.Path(cwd,directory))
            self.load_shifts(path)
            self.load_names(path)


    def dump(self):
        for key in sorted(self.shifts):
            print(key)
            print(self.names[key])
            print(self.shifts[key])
            print()

    def match(self, *target_shifts_sets):
        for i,target_shifts in enumerate(target_shifts_sets):

            shifts_str = ["%7.3f" % shift  for shift in target_shifts]
            print (f"set {i}: {', '.join(shifts_str)}")
            print()
            scores = {}
            # print(target_shifts)

            for id in self.shifts:
                # print(self.shifts[id])
                score =0.0
                for target_shift in target_shifts:
                    diffs  = [abs(shift - target_shift) for shift in self.shifts[id]]
                    score+=min(diffs)
                scores[score] = id

            sorted_scores = sorted(scores.keys())

            print(f"score    id             molecule")
            print(f"-----    --             --------")

            for i in range(10):
                best = scores[sorted_scores[i]]
                print("%-7.3f  %-10s     %-s" % (sorted_scores[i],best[0],self.names[best[0]].strip(' "')))




    def run(self, *shift_sets):
        self.load()
        for shift_set in shift_sets:
            self.match(shift_set)
        # self.dump()


if __name__ == '__main__':
    match = Match()

    shifts = [float(arg) for arg in sys.argv[1:]]

    match.run(shifts)