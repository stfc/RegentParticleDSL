import os
import subprocess
import sys
from datetime import datetime

x_counts = []
y_counts = []
z_counts = []
cutoffs = []

#Creates a series of tests with 1/3/5/7/9 particles in each dimension, and a variety of cutoff lengths
def define_tests():
    global x_counts, y_counts, z_counts, cutoffs
    x_counts=[1, 3, 6, 9]
    y_counts=[1, 3, 6, 9]
    z_counts=[1, 3, 6, 9]
#    for i in range(1,10,2):
#        x_counts.append(i)
#    for j in range(1,10,2):
#        y_counts.append(j)
#    for k in range(1,10,2):
#        z_counts.append(k)
    cutoffs.append(1.01)
#--    cutoffs.append(1.25)
#--    cutoffs.append(2)

#Runs the generator, solution computer and tests scripts on each of the inputs defined in define_tests
def run_tests():
    global x_counts, y_counts, z_counts, cutoffs
    define_tests()
    test_num = 0
    start=datetime.now()
    copy = subprocess.run(["cp", "examples/interaction_count/defaults.rg", "."])
    for i in range(len(x_counts)):
        for j in range(len(y_counts)):
            for k in range(len(z_counts)):
                for l in range(len(cutoffs)):
                    test_num = test_num + 1
                    print("--------------------------")
                    print("RUNNING TEST NUMBER {}".format(test_num))
                    print("--------------------------")
                    print(x_counts[i], y_counts[j], z_counts[k], cutoffs[l])
                    box_x = max(x_counts[i] * 1.0, cutoffs[l] * 2.0)
                    box_y = max(y_counts[j] * 1.0, cutoffs[l] * 2.0)
                    box_z = max(z_counts[k] * 1.0, cutoffs[l] * 2.0)
                    print(box_x, box_y, box_z)
                    try:
                        bop = subprocess.run(["python3", "tests/interaction_count/generate_regular.py", "--output", "tests/interaction_count/test_0.hdf5"
                                           , 
                                           "--x_count","{}".format(x_counts[i]), "--y_count", "{}".format(y_counts[j]),
                                           "--z_count", "{}".format(z_counts[k]), "--cutoff", "{}".format(cutoffs[l])], check=True,
                                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    except:
                        print("Failed for {}".format(bop.args))
                        sys.exit(1)
                    
                    try:
                        pop = subprocess.run(["python3", "tests/interaction_count/compute_solution.py", "--input", "tests/interaction_count/test_0.hdf5"
                                            , "--output", "tests/interaction_count/solution_0.hdf5",
                                            "--box_x", "{}".format(box_x),
                                            "--box_y", "{}".format(box_y), "--box_z", "{}".format(box_z)], check=True,
                                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

                    except:
                        print("Failed for {}".format(pop.args))
                        sys.exit(1)
#                    try:
                    pop3 = subprocess.run(["legion/language/regent.py", "tests/interaction_count/interaction_test_sym.rg", "-input", "tests/interaction_count/test_0.hdf5"
                                                                    , "-solution", "tests/interaction_count/solution_0.hdf5",
                                                        "-x_cell", "{}".format(box_x * 1.0),
                                                        "-y_cell", "{}".format(box_y), "-z_cell", "{}".format(box_z)], check=False,
                                                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    if pop3.returncode != 0:
                        print("Failed symmetric test {} {} {} {}".format(test_num, x_counts[i], y_counts[j], z_counts[k] ))
                        print("Sol args {}".format(pop.args))
                        print("Args {}".format(pop3.args))
                        sys.exit(1)
                    
                    pop4 = subprocess.run(["legion/language/regent.py", "tests/interaction_count/interaction_test_asym.rg", "-input", "tests/interaction_count/test_0.hdf5"
                        , "-solution", "tests/interaction_count/solution_0.hdf5",
                                                        "-x_cell", "{}".format(box_x * 1.0),
                                                        "-y_cell", "{}".format(box_y), "-z_cell", "{}".format(box_z)], check=False,
                                                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    if pop4.returncode != 0:
                        print("Failed asymmetric test {} {} {} {}".format(test_num,  x_counts[i], y_counts[j], z_counts[k]))
                        print("Sol args {}".format(pop.args))
                        print("Args {}".format(pop4.args))
                        sys.exit(1)
    print("Tests took {}".format( datetime.now()-start))



if __name__ == '__main__':
        run_tests()
