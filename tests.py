import os
import subprocess
import sys

x_counts = []
y_counts = []
z_counts = []
cutoffs = []

def define_tests():
    global x_counts, y_counts, z_counts, cutoffs
    for i in range(1,10):
        x_counts.append(i)
    for j in range(1,10):
        y_counts.append(j)
    for k in range(1,10):
        z_counts.append(k)
    cutoffs.append(1.75)


def run_tests():
    global x_counts, y_counts, z_counts, cutoffs
    define_tests()
    test_num = 0
    for i in range(len(x_counts)):
        for j in range(len(y_counts)):
            for k in range(len(z_counts)):
                for l in range(len(cutoffs)):
                    print(x_counts[i], y_counts[j], z_counts[k], cutoffs[l])
                    box_x = max(x_counts[i] * 1.0, cutoffs[l] * 2.0)
                    box_y = max(y_counts[j] * 1.0, cutoffs[l] * 2.0)
                    box_z = max(z_counts[k] * 1.0, cutoffs[l] * 2.0)
                    print(box_x, box_y, box_z)
                    try:
                        bop = subprocess.run(["python3", "tests/interaction_count/generate_regular.py", "--output", "tests/interaction_count/test_{}.hdf5"
                                           .format(test_num), 
                                           "--x_count","{}".format(x_counts[i]), "--y_count", "{}".format(y_counts[j]),
                                           "--z_count", "{}".format(z_counts[k]), "--cutoff", "{}".format(cutoffs[l])], check=True)
                    except:
                        print("Failed for {}".format(bop.args))
                        sys.exit(1)
                    
                    try:
                        pop = subprocess.run(["python3", "tests/interaction_count/compute_solution.py", "--input", "tests/interaction_count/test_{}.hdf5"
                                            .format(test_num), "--output", "tests/interaction_count/solution_{}.hdf5".format(test_num),
                                            "--box_x", "{}".format(box_x),
                                            "--box_y", "{}".format(box_y), "--box_z", "{}".format(box_z)], check=True)

                    except:
                        print("Failed for {}".format(pop.args))
                        sys.exit(1)
#                    try:
                    pop3 = subprocess.run(["regent", "tests/interaction_count/interaction_test_asym.rg", "-input", "tests/interaction_count/test_{}.hdf5"
                                                                    .format(test_num), "-solution", "tests/interaction_count/solution_{}.hdf5".format(test_num),
                                                        "-x_cell", "{}".format(box_x * 1.0),
                                                        "-y_cell", "{}".format(box_y), "-z_cell", "{}".format(box_z)], check=False)
                    if pop3.returncode != 0:
                        print("Failed test {}".format(test_num))
                        sys.exit(1)



if __name__ == '__main__':
        run_tests()
