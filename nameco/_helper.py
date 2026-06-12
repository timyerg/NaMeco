import os
import subprocess
from importlib.metadata import version


#Function to run bash commands
def bash(cmd):
    return subprocess.check_output(cmd, shell=True)


#Function to print greetings 
def greetings(log, LOGS):
    bash(f'mkdir -p {LOGS}')
    tool_version = version("nameco")
    grt = """
######################################################
#                                                    #
# ##     #         ##     ##                         #                
# # #    #         # #   # #                         #
# # #    #         #  # #  #                         #
# #  #   #  #####  #   #   #   ###     ###     ###   #
# #  #   # #     # #       #  #   #   #   #   #   #  #
# #   #  #       # #       # #     # #       #     # #
# #   #  #  ###### #       # ####### #       #     # #
# #    # # #     # #       # #       #       #     # #
# #    # # #     # #       #  #   #   #   #   #   #  #
# #     ##  #####  #       #   ###     ###     ###   #
#                                                    #
######################################################

Powered by Coffee
You are running NaMeco v{}
If you used this pipeline, please cite our paper: https://doi.org/10.1186/s12864-025-12415-x
Also, don't forget to cite the tools that were used in this pipeline
    """
    print(grt.format(tool_version))
    bash(f'echo "You are running NaMeco v{tool_version}" >> {log}')
    

#Function to wrap messages into hashtags
def hashtags_wrapper(sub):
    print(f"\n{'#'*(len(sub)+8)}\n### {sub} ###\n{'#'*(len(sub)+8)}\n")
    
    
#Function to check logs
def log_checker(log, samples, file, pat='done. Enjoy\n'):
    skip, checks = [], []
    if os.path.exists(log):
        with open(log, 'rt') as txt:
            skip = [l.split(' ')[0] for l in txt.readlines() if l.endswith(pat)]
    for sample in samples:
        checks.append(os.path.exists(file.format(sample)))
        checks.append(sample in skip)
    return skip, checks
