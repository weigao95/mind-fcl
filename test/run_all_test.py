import os
import sys
import multiprocessing
import subprocess


def print_with_color(string, color, highlight=False):
    """
    Colored print
    colorlist:
        red,green
    """
    end = "\033[1;m"
    pstr = ""
    if color == "red":
        if highlight:
            pstr += '\033[1;41m'
        else:
            pstr += '\033[1;31m'
    elif color == "green":
        if highlight:
            pstr += '\033[1;42m'
        else:
            pstr += '\033[1;32m'
    elif color == "yellow":
        if highlight:
            pstr += '\033[1;43m'
        else:
            pstr += '\033[1;33m'

    else:
        print(("Error Unsupported color:" + color))

    if isinstance(string, str):
        print((pstr + string + end))
    else:
        print((pstr + string.decode('utf-8') + end))


def run_single_test(test_binary: str):
    """
    Run a single test, log the info if failed
    :param test_binary:
    :return:
    """
    try:
        print('Trying {test_case}'.format(test_case=test_binary))
        output = subprocess.check_output(test_binary, shell=False)
        return True, output
    except subprocess.CalledProcessError as e:
        output = e.output
        return False, output


def guess_bin_directory(release_only: bool):
    """
    Guess a list of binary directories that may contains the tests
    :return:
    """
    test_dir = os.path.dirname(os.path.realpath(__file__))
    project_root = os.path.join(test_dir, os.pardir)

    # Try with build/cmake-build-debug/cmake-build-release
    build_bin_dirs = list()
    attempt_dir_list = ['build', 'cmake-build-release', 'out/build/x64-Release']
    if not release_only:
        attempt_dir_list.append('cmake-build-debug')

    # Find the binary
    for dir_attempt in attempt_dir_list:
        build_dir = os.path.join(project_root, dir_attempt)
        build_dir = os.path.abspath(build_dir)
        build_bin_dir = os.path.join(build_dir, 'test')
        if os.path.exists(build_bin_dir):
            build_bin_dirs.append(os.path.abspath(build_bin_dir))

    # OK
    return build_bin_dirs


def collect_test_binaries(bin_dir_list):
    """
    Given the list of binary directories, collect the test cases
    :param bin_dir_list:
    :return:
    """
    binary_path_list = list()
    for binary_dir in bin_dir_list:
        dir_contents = os.listdir(binary_dir)
        for filename in dir_contents:
            is_binary = False
            if sys.platform == "win32":
                is_binary = filename.startswith('test_') and filename.endswith('.exe')
            else:  # Linux or mac
                is_binary = filename.startswith('test_') and (not filename.endswith('.cmake'))

            # App to path list if this is a binary
            if is_binary:
                binary_path = os.path.join(binary_dir, filename)
                binary_path_list.append(str(binary_path))
    return binary_path_list


def run_tests(release_only: bool = True, parallel: bool = False):
    # Get the test cases
    bin_dirs = guess_bin_directory(release_only)
    binary_test_cases = collect_test_binaries(bin_dirs)
    if len(binary_test_cases) == 0:
        print_with_color('Cannot find any test cases, thus exit.', "red")
        return

    # Try pool based solution
    if parallel:
        n_process = 32
        pool = multiprocessing.Pool(processes=n_process)
        output_list = pool.map(run_single_test, binary_test_cases)
    else:
        output_list = list()
        for i in range(len(binary_test_cases)):
            output_i = run_single_test(binary_test_cases[i])
            output_list.append(output_i)

    # Output the result
    for i in range(len(output_list)):
        passed_i, output_i = output_list[i]
        if not passed_i:
            print_with_color(output_i, 'red')


if __name__ == '__main__':
    run_tests(release_only=True)
