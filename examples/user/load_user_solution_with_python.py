import json, numpy as np

def load_psopt_solution(path):
    sol = json.load(open(path))
    for p in sol["phases"]:
        for k in ("time", "states", "controls", "costates", "relative_errors", "parameters"):
            p[k] = np.asarray(p[k])     # states/controls/costates -> (n_vars, n_nodes)
    return sol

sol = load_psopt_solution("user_solution.json")
ph1 = sol["phases"][0]
t, x, u, lam = ph1["time"], ph1["states"], ph1["controls"], ph1["costates"]
# x[i] is the trajectory of state i+1; u[j] of control j+1
