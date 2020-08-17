import networkx as nx
import pickle
import os
import copy
import argparse
from collections import namedtuple

from parsers import parse_grid_types, parse_pattern,\
                    parse_vpr_delay_table, parse_pattern_delays
from feeder import TimingGraph, Stretcher, Timer
from ilp_placer import ILP
from ilp_placer_types import *
from selector import Selector
from legalizer import Legalizer, Lut


SCALEUP = 1e12
#Multiplicative constant used to scale delays before integer casting
W_sel = 1
#Movement region seen by the selection LP
W_stretch = 1
#Movement region seen by the placement ILP 
TMLIM = 600
#Solver time-out in seconds
MIN_IMP = 10
#Improvement threshold used on edges when selecting movable nodes

CHECKS_ON = True
#Enables timing consistency checks
TOLERANCE = 10
#Tolerance for the consistency mismatch.
#NOTE: Must be increased for deep benchmarks such as LU8PEEng,
#as the cumulative rounding error will be larger than 10 ps.

ILP_RELAXATION = 1.1
#Maximum relaxation of the ILP target delay, with respect to the LP target delay
BIN_SEARCH_STOP = 0.005
#Controls when the binary search for the target delay stops
PASS_LB = False
#If true, the critical path delay is bounded from below based on the previous failures
STORE_STATE = True
#If true, stores the state after each major step in the flow 

parser = argparse.ArgumentParser()
parser.add_argument("benchmark")
parser.add_argument("seed")
args = parser.parse_args()

benchmark = args.benchmark
seed = int(args.seed)

import_dir = "import_dir/"
import_subdir = import_dir + "%s/%d/" % (benchmark, seed)
ftg_file = import_subdir + "ftg.pickle"
grid_type_file = import_subdir + "grid_types.dump"
pat_file = import_dir + "best_pattern_i13.txt"
grid_delay_file = import_subdir + "lookup_dump.echo"
pat_delay_file = import_dir + "delays.dump"

grid_types = parse_grid_types(grid_type_file)
grid_delays = parse_vpr_delay_table(grid_delay_file)
pattern = parse_pattern(pat_file)
pattern_delays = parse_pattern_delays(pattern, pat_delay_file)

with open(ftg_file, "r") as dumpf:
    ftg = pickle.load(dumpf)

ftg_upscaled = copy.deepcopy(ftg)
ftg_upscaled.scale_delays(SCALEUP, ftg_upscaled.graph)
ftg_upscaled.sta(ftg_upscaled.graph)
sel_timer = Timer(ftg.graph, pattern_delays, grid_delays, grid_types)
sel_timer.set_scaleup(SCALEUP)
sel = Selector(ftg_upscaled.graph, sel_timer, W_sel)

State = namedtuple("State", ["core", "coverage", "legalized_tg",\
                             "orig_cpd", "cov_cpd", "leg_cpd"])

##########################################################################
def iter_ilp(sel, ftg, target_cpd, fix_covs = [], cpd_lb = None,\
             incumbent_coverage = []):
    """Performs one full ILP construction and solving iteration.

    Parameters
    ----------
    sel : Selector
        The selector object holding the core.
    ftg : feeder.TimingGraph
        The feeder timing graph holding the compressed timing graph,
        around the edges of the core. This should be the version with
        appropriately scaled delays.
    target_cpd : float
        The target critical path delay.
    fix_covs : Optional[Dict or List]
        Specifies coverage to be enforced in the solution.
    forced_cpd_lb : Optional[float]
        Optional CPD lower bound. This is useful when solving
        a sequence of related ILPs (as is the case here with the binary search).
    incumbent_coverage : Optional[List[Tuple[str]]]
        Optional incumbent solution coverage. Needed if coverage similarity is
        to be maximized.
    
    Returns
    -------
    Tuple[Dict[Tuple[str], PhysDirCon], float] or Tuple[None, float]
        If a feasible solution is found, returns a dictionary of covered
        edges with their covering pairs, along with the achieved CPD.
        Otherwise, returns the new lower bound on CPD
    """
 
    #Find maximum edge delays based on slack
    #(If relaxation is reasonably high, does not prune much)
    sel.cpd = target_cpd
    max_tds = {}
    max_tds_str = {}
    max_sel_tds = sel.sta_bound_delays()
    for e in max_sel_tds:
        u, v = e
        try:
            uc = ftg.node_counts[u]
        except:
            for p in ftg.compressed_graph.pred[v]:
                if p.startswith("dummy_"):
                    orig = ftg.compressed_graph.node[p].get("orig", None)
                    if orig == u:
                        uc = ftg.node_counts[p]
                        u = p
                        break
        try:
            vc = ftg.node_counts[v]
        except:
            for c in ftg.compressed_graph[u]:
                if c.startswith("dummy_"):
                    orig = ftg.compressed_graph.node[c].get("orig", None)
                    if orig == v:
                        vc = ftg.node_counts[c]
                        v = c
                        break
        eb = BasicEdge(uc, vc)
        max_tds.update({eb : max_sel_tds[e]})
        max_tds_str.update({(u, v) : max_sel_tds[e]})

    #Export the cov and stretch maps.
    stc.max_tds = max_tds_str
    cov_map = stc.export_cov_map(ftg.node_counts)
    stretch_map = stc.export_stretch_map(ftg.node_counts)

    #Prune edges that can never become critical.
    max_achievable_delays = {}
    min_achievable_delays = {}
    inv_node_counts = {ftg.node_counts[u] : u[:-1] + 'd'\
                       if u.endswith(".q") else u\
                       for u in ftg.node_counts}

    cov_map, stretch_map = stc.resolve_node_conflicts(cov_map, stretch_map, inv_node_counts)

    for e in set(cov_map.keys() + stretch_map.keys()):
        #NOTE: For this to work, the feeder must not filter the slow pairs.
        #This must be left to ilp_placer. 
        u, v = e
        u_str = inv_node_counts[u]
        if u_str.endswith(".d"):
            u_str = u_str[:-1] + 'q'
        if not ftg.graph.has_node(u_str):
            u_str = ftg.compressed_graph.node[u_str]["orig"]
        v_str = inv_node_counts[v]
        if v_str.endswith(".q"):
            v_str = v_str[:-1] + 'd'
        if not ftg.graph.has_node(v_str):
            v_str = ftg.compressed_graph.node[v_str]["orig"]
        maad = max([c.td for c in cov_map.get(e, [])]\
            + [s.td for s in stretch_map.get(e, [])])
        max_achievable_delays.update({(u_str, v_str) : maad})
        miad = min([c.td for c in cov_map.get(e, [])]\
            + [s.td for s in stretch_map.get(e, [])])
        min_achievable_delays.update({(u_str, v_str) : miad})
    
    max_achievable_delays = stc.export_max_achievable_delays()
    irrelevant = sel.prune_irrelevant_edges(max_achievable_delays)
    print("irrelevant: %d" % len(irrelevant))
    for u, v in irrelevant:
        u_int = ftg.node_counts[u]
        v_int = ftg.node_counts[v]
        e = BasicEdge(u_int, v_int)
        try:
            del cov_map[e]
        except:
            pass
        try:
            del stretch_map[e]
        except:
            pass 

    #Prepare inputs for ILP construction.
    fixed_nodes = {ftg.node_counts[u]\
                  : ftg.compressed_graph.node[u]["coords"]\
                    for u in ftg.fixed}
    
    init_positions = {ftg.node_counts[u]\
                     : ftg.compressed_graph.node[u]["coords"]\
                       for u in ftg.compressed_graph}

    #Construct the ILP.
    ilp = ILP(ftg.export_tg(), ftg.export_tarrs(),\
              ftg.export_node_counts(), init_positions,\
              fixed_nodes, cov_map, stretch_map, max_tds,\
              max_achievable_delays, sel.cpd, cpd_lb)
    ilp.load_pattern_edges([(e.u, e.v) for e in pattern])
    ilp.build_all()
    aliases = {ftg.node_counts[u] : u for u in ftg.node_counts}
    ilp.alias_nodes(aliases)
    ilp.generate_constraints()
    ilp.set_incumbent_coverage(incumbent_coverage)

    #Solve the ILP.
    ilp.solve(TMLIM)
    
    #Parse the solution.
    sol = ilp.parse_solution()
    if not sol:
        return [None, target_cpd]#ilp.get_lb()]

    ilp.get_arrival_times()
    ilp.get_edge_delays(TOLERANCE)
    ilp.get_positions()
    if CHECKS_ON:
        ilp.check_solution_timing_consistency(TOLERANCE)
        ilp.check_cov_consistency(pattern)
    coverage = ilp.reconstruct_coverage()

    return [coverage, ilp.solution["Tamax"]]
##########################################################################

sel.find_bound()
print "Min CPD: ", sel.cpd
absolute_min_outer_cpd = sel.cpd
min_outer_cpd = sel.cpd
initial_cpd = ftg.sta(ftg.graph)
print "Initial CPD:", initial_cpd
max_outer_cpd = initial_cpd * SCALEUP
absolute_best = max_outer_cpd
absolute_best_core = None
absolute_best_placement = {}

freeze = False
while max_outer_cpd - min_outer_cpd > max_outer_cpd * BIN_SEARCH_STOP:
    target_outer_cpd = 0.5 * (min_outer_cpd + max_outer_cpd)
    relaxation = float(target_outer_cpd) / absolute_min_outer_cpd
    print target_outer_cpd

    #Find a minimum improvement set that meets the relaxed bound.
    sel.set_relaxation(relaxation)
    sel.select_edges_for_improvement(MIN_IMP, TMLIM)
    if not len(sel.min_edge_improvements):
        print("LP infeasible")
        min_outer_cpd = target_outer_cpd
        continue
    print("cpd = %.4g, |V| = %d, |E| = %d"\
          % (sel.cpd / SCALEUP, len(sel.affected_nodes),\
             len(sel.min_edge_improvements)))

    #Compress the feeder graph.
    ftg.find_affected_edges(sel.min_edge_improvements.keys())
    ftg.compress_tg()
    ftg.count_nodes()
    ftg.push_memory_dummies()
 
    #Prepare the Stretcher inputs.
    ftg_upscaled = copy.deepcopy(ftg)
    ftg_upscaled.scale_delays(SCALEUP, ftg_upscaled.graph)
    ftg_upscaled.sta(ftg_upscaled.graph)
    ftg_upscaled.scale_delays(SCALEUP, ftg_upscaled.compressed_graph)
    ftg_upscaled.sta(ftg_upscaled.compressed_graph)
    timer = Timer(ftg.compressed_graph, pattern_delays, grid_delays, grid_types)
    timer.set_scaleup(SCALEUP)
    init_positions = {u : ftg_upscaled.compressed_graph.node[u]["coords"]\
                          for u in ftg_upscaled.compressed_graph}
    
    #Initialize the Stretcher.
    stc = Stretcher(ftg_upscaled.cov, ftg_upscaled.stretch, ftg_upscaled.fixed,\
                    init_positions, pattern, grid_types, W_stretch, {}, timer)
    stc.cluster_dir_cons()
 
    #Inner loop: iteratively relaxes the ILP bound until feasibility
    #is reached or the target delay reaches maximum tolerance before
    #we seek a new core.

    orig_sel_cpd = sel.cpd
    target_cpd = sel.cpd
    min_target_cpd = sel.cpd
    max_target_cpd = min(sel.cpd * ILP_RELAXATION, absolute_best)
   
    found_witness = False
    while max_target_cpd - min_target_cpd > max_target_cpd * BIN_SEARCH_STOP:
        target_cpd = 0.5 * (min_target_cpd + max_target_cpd)
        print "target:", target_cpd
        coverage, cpd = iter_ilp(sel, ftg_upscaled, target_cpd,\
                                 cpd_lb = (min_target_cpd if PASS_LB else None),\
                                 incumbent_coverage = list(absolute_best_placement.keys()))
        if coverage is None:
            min_target_cpd = cpd
            print "failed:", cpd
        else:
            found_witness = True
            #NOTE: If irrelevant edge pruning is enabled, we must not rely on better-than
            #target delays, as they may not be achievable with all edges present.
            cpd = target_cpd
            max_target_cpd = cpd
            print "achieved:", cpd
            if cpd < absolute_best:
                absolute_best = cpd
                print "incumbent:", absolute_best
                absolute_best_core = copy.deepcopy(sel)
                print "Previously covered:", len(absolute_best_placement)
                print "New covered:", len(coverage)
                print "Intersection:",\
                len(set(coverage.keys()).intersection(set(absolute_best_placement.keys())))
                absolute_best_placement = coverage
        os.system("rm -f ilp.sol")
        #if found_witness and not freeze:
        #    break
    sel.cpd = orig_sel_cpd
    if found_witness:
        max_outer_cpd = max_target_cpd
    else:
        min_outer_cpd = min_target_cpd
    
print absolute_best

#NOTE: Keep the intermediate state in case legalization fails.
if STORE_STATE:
    state = State(core = absolute_best_core, coverage = absolute_best_placement,\
                  legalized_tg = None, orig_cpd = initial_cpd,\
                  cov_cpd = -1, leg_cpd = -1)
    with open("state_%s_%d_Wsel_%d_Wstretch_%d.dump"\
              % (benchmark, seed, W_sel, W_stretch), "w") as dumpf:
        pickle.dump(state, dumpf)

with open("state_%s_%d_Wsel_%d_Wstretch_%d.dump"\
          % (benchmark, seed, W_sel, W_stretch), "r") as dumpf:
    state = pickle.load(dumpf)

absolute_best_core = state.core
absolute_best_placement = state.coverage

if not absolute_best_placement:
    print "Nothing to legalize."
    exit(0)

legalizer = Legalizer(ftg.graph, grid_types, absolute_best_placement,\
                      pattern_delays, grid_delays)
legalizer.refresh_timing()
cov_cpd = legalizer.sta()
print "CPD after coverage:", cov_cpd
legalizer.compute_criticalities()
legalizer.build_clusters()
print "Overflowed clusters:"
for cls in legalizer.clusters:
    ovf = legalizer.clusters[cls].overflowed()
    if ovf > 0:
        print cls, ovf
legalizer.legalize()
legalizer.commit_to_tg()
print "Post-legalization:"
print "Overflowed clusters:"
for cls in legalizer.clusters:
    ovf = legalizer.clusters[cls].overflowed()
    if ovf > 0:
        print cls, ovf
leg_cpd = legalizer.sta()
print "CPD after legalization:", leg_cpd

if STORE_STATE:
    state = State(core = absolute_best_core, coverage = absolute_best_placement,\
                  legalized_tg = legalizer.tg, orig_cpd = initial_cpd,\
                  cov_cpd = cov_cpd, leg_cpd = leg_cpd)
    with open("state_%s_%d_Wsel_%d_Wstretch_%d.dump"\
              % (benchmark, seed, W_sel, W_stretch), "w") as dumpf:
        pickle.dump(state, dumpf)

with open("cls.out", "w") as outf:
    for i, cls in enumerate(sorted(legalizer.clusters)):
        txt = "Cluster %d (%d, %d)\n" % (i, cls[0], cls[1])
        line = '-' * len(txt) + "\n\n"
        txt += line
        for lut in sorted(legalizer.clusters[cls].luts, key = lambda l : (l.pos, l.o)):
            txt += "%d %s\n" % (lut.pos, lut.o)
        txt += line
        outf.write(txt)
