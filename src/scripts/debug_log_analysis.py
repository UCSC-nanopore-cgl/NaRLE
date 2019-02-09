from __future__ import print_function
import argparse
import glob
import numpy as np

class debug_attributes:
    def __init__(self, i, line_parts):
        self.function = ''
        self.uuid = None
        self.debug_type = ''
        self.value = None
        self.extra_info = None

        parts_len = len(line_parts)
        self.debug_type = line_parts[i].strip()
        if i-1 >= 0: self.function = line_parts[i-1].strip()
        if i-2 >= 0: self.uuid = line_parts[i-2].strip()
        if i+1 < parts_len: self.value = line_parts[i+1].strip()
        if i+2 < parts_len: self.extra_info = line_parts[i+2:]

def main(args):

    input_logs = glob.glob(args.input)
    if len(input_logs) == 0:
        print("Could not find files matching {}".format(args.input))
        return

    all_log_data = dict()
    all_log_data_by_type = dict()
    all_log_data_by_function = dict()

    for input_log in input_logs:
        # uuid -> function -> debug key -> value -> extra values
        log_data = list()
        log_data_by_type = dict()
        log_data_by_function = dict()
        all_log_data[input_log] = log_data
        all_log_data_by_type[input_log] = log_data_by_type
        all_log_data_by_function[input_log] = log_data_by_function

        def save_attrs(log_attrs):
            if log_attrs.function not in log_data_by_function: log_data_by_function[log_attrs.function] = list()
            if log_attrs.debug_type not in log_data_by_type: log_data_by_type[log_attrs.debug_type] = list()

            log_data.append(log_attrs)
            log_data_by_type[log_attrs.debug_type].append(log_attrs)
            log_data_by_function[log_attrs.function].append(log_attrs)

        # read in
        with open(input_log, 'r') as log:
            last_line = None
            for raw_line in log:
                # only get debug lines, only get them once
                if "DEBUG_" not in raw_line: continue
                tmp_line = last_line
                last_line = raw_line
                if tmp_line is not None and (raw_line.endswith(tmp_line) or tmp_line.endswith(raw_line)): continue

                # something is debugged in here
                line_parts = list(map(lambda x: x.strip(), raw_line.split(":")))
                for i, part in enumerate(line_parts):
                    if part.startswith("DEBUG_"):
                        attrs = debug_attributes(i, line_parts)
                        save_attrs(attrs)
                        break


        # functions and stats
        functions = list(log_data_by_function.keys())
        types = list(log_data_by_type.keys())
        functions.sort()
        types.sort()

        # gather stats
        all_t = dict()
        all_f = dict()
        f_to_t = dict()
        t_to_f = dict()
        for d in log_data:
            t = d.debug_type
            f = d.function
            v = d.value

            if t not in all_t:
                all_t[t] = list()
                t_to_f[t] = dict()
            if f not in all_f:
                all_f[f] = list()
                f_to_t[f] = dict()
            if f not in t_to_f[t]: t_to_f[t][f] = list()
            if t not in f_to_t[f]: f_to_t[f][t] = list()

            all_t[t].append(v)
            all_f[f].append(v)
            t_to_f[t][f].append(v)
            f_to_t[f][t].append(v)

        # print stats
        def print_stats(values, tabs=0):
            try:
                int_values = list(map(int, values))
                print("{}MAX  : {}".format("\t" * tabs, max(int_values)))
                print("{}MIN  : {}".format("\t" * tabs, max(int_values)))
                print("{}MEAN : {}".format("\t" * tabs, np.mean(int_values)))
                print("{}STD  : {}".format("\t" * tabs, np.std(int_values)))
            except:
                print("{}Do not appear to be all integers".format("\t" * tabs))

        print("\nBy Type:")
        for t in types:
            print("\t{}:".format(t))
            for f in functions:
                if f in t_to_f[t]:
                    print("\t\t{}:".format(f))
                    print_stats(t_to_f[t][f], 3)
        print("\nBy Function:")
        for f in functions:
            print("\t{}:".format(f))
            for t in types:
                if t in f_to_t[f]:
                    print("\t\t{}:".format(t))
                    print_stats(f_to_t[f][t], 3)



        # try:
        #     if v is None or int(v) < 0: continue
        # except:  # v is not of type int
        #     continue
        # v = int(v)







if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--input', '-i', dest='input', required=True, type=str, help="input file (glob)")

    args = parser.parse_args()
    main(args)