import os
from pandas import cut
import shared
# Assume that the codes are created

runMatlabCmd = """
matlab -nodisplay -nojvm -nosplash -nodesktop -sd "/home/lev/code/research/pfm/matlab/" -r "try, run('/home/lev/code/research/pfm/matlab/get_block_error.m'), catch, exit(1), end, exit(0);"
echo "matlab exit code: $?"
"""

error_ps = [0.005, 0.01, 0.02, 0.03, 0.04]
ks = range(10, 14)
Ntrials = 100
cutoff = 100
for p in error_ps:
    regularity = (5, 6)
    for i in ks:
        graph_file_name = shared.get_graph_file_name(i, regularity)
        paramsf = open('../matlab/tmp/params.json', "w")
        paramsf.write(
            "{" +
            f'"p": {p}, "n_trials": {Ntrials}, "cutoff": {cutoff}, "file_out": "../matlab/out/{i * regularity[1]}-{i * regularity[0]}-graph-{p}-{Ntrials}-{cutoff}.json"' + "}"
        )
        os.system(runMatlabCmd)
				# print(runMatlabCmd
