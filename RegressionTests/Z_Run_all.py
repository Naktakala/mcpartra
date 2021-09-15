import subprocess
import os
import sys

# To run a range of tests pass a text string via the command line.
# i.e., "tests_to_run=[*range(14,16)]"

# This python script executes the regression test suite.
# In order to add your own test, copy one of the test blocks
# and modify the logic to what you would like tested.

# General guidance: Do not write a regression test that checks
# for number of iterations. Rather make it check for answers
# since something as simple as a preconditioner can change
# iteration count but still produce the same answer

kscript_path = os.path.dirname(os.path.abspath(__file__))
kchi_src_pth = kscript_path + '/../'
kpath_to_exe = kchi_src_pth + '/bin/mcpartra'
tests_to_run = []
print_only   = False

print("")
if len(sys.argv) >= 2:
    exec(sys.argv[1])
print("************* MCParTra Regression Tests *************")
print("")
test_number = 0
num_failed = 0

# Determine if we are on TACC (each test will require a separate job)
hostname = subprocess.check_output(['hostname']).decode('utf-8')
tacc = True if "tacc.utexas.edu" in hostname else False

def format3(number):
    return "{:3d}".format(number)


def format_filename(filename):
    return "{:35s}".format(filename)

def parse_output(out, search_strings_vals_tols):
    global num_failed
    test_passed = True
    for search in search_strings_vals_tols:
        find_str = search[0]
        true_val = search[1]
        tolerance = search[2]

        # start of the string to find (<0 if not found)
        test_str_start = out.find(find_str)
        # end of the string to find
        test_str_end = test_str_start + len(find_str)
        # end of the line at which string was found
        test_str_line_end = out.find("\n", test_str_start)

        test_passed = True
        if test_str_start >= 0:
            # convert value to number
            test_val = float(out[test_str_end:test_str_line_end])
            if not abs(test_val - true_val) < tolerance:
                test_passed = False
        else:
            test_passed = False

    if test_passed:
        print(" - Passed")
    else:
        print(" - FAILED!")
        num_failed += 1
        print(out)

    return test_passed

def run_test_tacc(file_name, comment, num_procs,
		search_strings_vals_tols):
    test_name = format_filename(file_name) + " " + comment + " " + str(num_procs) + " MPI Processes"
    print("Running Test " + format3(test_number) + " " + test_name, end='', flush=True)
    if print_only: print(""); return
    with open(f"RegressionTests/{file_name}.job", 'w') as job_file:
        job_file.write(
f"""#!/usr/bin/bash
#
#SBATCH -J {file_name} # Job name
#SBATCH -o RegressionTests/{file_name}.o # output file
#SBATCH -e RegressionTests/{file_name}.e # error file
#SBATCH -p skx-normal # Queue (partition) name
#SBATCH -N {num_procs // 48 + 1} # Total # of nodes
#SBATCH -n {num_procs} # Total # of mpi tasks
#SBATCH -t 00:05:00 # Runtime (hh:mm:ss)
#SBATCH -A Massively-Parallel-R # Allocation name (req'd if you have more than 1)

ibrun {kpath_to_exe} RegressionTests/{file_name}.lua master_export=false"""
        )
    os.system(f"sbatch -W RegressionTests/{file_name}.job > /dev/null")  # -W means wait for job to finish
    with open(f"RegressionTests/{file_name}.o", 'r') as outfile:
        out = outfile.read()

    passed = parse_output(out, search_strings_vals_tols)

    # Cleanup
    if passed:
        os.system(f"rm RegressionTests/{file_name}.job RegressionTests/{file_name}.o RegressionTests/{file_name}.e")


def run_test_local(file_name, comment, num_procs,
		search_strings_vals_tols):
    test_name = format_filename(file_name) + " " + comment + " " + str(num_procs) + " MPI Processes"
    print("Running Test " + format3(test_number) + " " + test_name, end='', flush=True)
    if print_only: print(""); return
    process = subprocess.Popen(["mpiexec", "-np", str(num_procs), kpath_to_exe,
                                "RegressionTests/" + file_name + ".lua", "master_export=false"],
                               cwd=kchi_src_pth,
                               stdout=subprocess.PIPE,
                               universal_newlines=True)
    process.wait()
    out, err = process.communicate()

    parse_output(out, search_strings_vals_tols)

def run_test(file_name, comment, num_procs,
		search_strings_vals_tols):
    global test_number
    test_number += 1
    if ((tests_to_run) and (test_number in tests_to_run)) or \
       (not tests_to_run):
        if tacc:
            run_test_tacc(file_name, comment, num_procs,
                search_strings_vals_tols)
        else:
            run_test_local(file_name, comment, num_procs,
                search_strings_vals_tols)

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ MCParTra tests

run_test(
    file_name="MonteCarlo1D_1MatSource",
    comment="1D Test - Material source",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value1=", 2.02976 , 1.0e-5],
                              ["[0]  Max-value2=", 0.994566, 1.0e-5]])

run_test(
    file_name="MonteCarlo2D_1Poly",
    comment="2D Test - Material source",
    num_procs=1,
    search_strings_vals_tols=[["[0]  Max-value1=", 3.33722 , 1.0e-5],
                              ["[0]  Max-value2=", 0.249672, 1.0e-6]])


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END OF TESTS
print("")
if num_failed == 0:
    print("All regression tests passed!")
else:
    print("ERROR: Not all regression tests passed!")
    print("Number of Tests failed = " + str(num_failed))
print("")
print("************* End of Regression Test *************")
print("")
if num_failed == 0:
    sys.exit(0);
else:
    sys.exit(1); 
