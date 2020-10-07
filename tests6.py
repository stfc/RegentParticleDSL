import os
import subprocess
import sys
from datetime import datetime

#Runs the generator, solution computer and tests scripts on each of the inputs defined in define_tests
def run_tests():
    global x_counts, y_counts, z_counts, cutoffs
    define_tests()
    test_num = 0
    start=datetime.now()
    copy = subprocess.run(["cp", "examples/interaction_count/defaults.rg", "."])
    pop3 = subprocess.run(["legion/language/regent.py", "tests/zero_part/zero_part.rg"], check=False,
                                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if pop3.returncode != 0:
        print("Failed zero part test")
        print("Args {}".format(pop3.args))
        sys.exit(1)
    
    print("Test took {}".format( datetime.now()-start))



if __name__ == '__main__':
        run_tests()
