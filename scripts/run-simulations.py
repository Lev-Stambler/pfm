from pandas import cut
import shared
# Assume that the codes are created

error_ps = [0.005, 0.01, 0.02, 0.03, 0.04]
ks = range(5, 8)
Ntrials = 1000
cutoff = 30
for i in ks:
    regularity = (5, 6)
    graph_file_name = shared.get_graph_file_name(i, regularity)
    for p in error_ps:
        paramsf = open('../matlab/tmp/params.json', "w")
        paramsf.write(
            "{" + f'"p": {p}, "n_trials": {Ntrials}, "cutoff": {cutoff}, "file_out": "../matlab/out/{i * regularity[1]}-{i * regularity[0]}-graph-{p}-{Ntrials}-{cutoff}.json"' + "}"
        )
