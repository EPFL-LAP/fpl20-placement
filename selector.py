"""Module contains routines for selecting the edges for improvement.
"""

import networkx as nx
import os
import copy
import time

from ilp_placer import TimingGraph
from ilp_placer_types import check_nonneg_float, check_nonneg_int, NodePos, NodeCls,\
                             PhysDirCon, PhysProgCon, FreeDirCon, FreeProgCon
from feeder import Timer

precision_4 = lambda f : "%.4f" % f
precision_0 = lambda f : "%d" % f
print_float = precision_0

OPTIMIZER = "cplex"

##########################################################################
class ImprovementGraph(TimingGraph):
    """Extends the TimingGraph class of the ilp_placer module so
    that the edge delays are represented by a fixed value and a variable
    improvement, which models the edge-selection process.

    Paramters
    ---------
    *TimingGraph

    Attributes
    ----------
    *TimingGraph

    Methods
    -------
    *TimingGraph
    generate_td_imp_constraints()
        Generates the improvement constraints.
    generate_td_imp_bounds(imp_bounds : Dict[Tuple[str], float],
                           positive : bool )
        Generates the improvement bounds.
    generate_constraints()
        Redefined to include new constraints.
    print_constraints()
        Redefined to print new constraints.
    print_bounds()
        Redefined to print new bounds.
    """

    #-----------------------------------------------------------------------#
    def __init__(self, tg, ch, tarrs, node_counts):
        """Constructor of the ImprovementGraph class.
        """

        super(ImprovementGraph, self).__init__(tg, ch, tarrs, node_counts)
        self.imp_bounds = {}

        self.edge_delay_constraints = []
        self.edge_delay_bounds = []
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def generate_td_imp_constraints(self):
        """Generates the improvement constraints.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        #String conversions:
        delay_var = lambda u, v : "td_%d_%d" % (self.node_counts[u],\
                                                self.node_counts[v])
        imp_var = lambda u, v : "imp_%d_%d" % (self.node_counts[u],\
                                               self.node_counts[v])

        for u, v in self.ch:
            cst = delay_var(u, v) + " + " + imp_var(u, v) + " >= "\
                + print_float(self.tg[u][v]["td"])
            self.edge_delay_constraints.append(cst)
            if self.imp_bounds:
                bound = imp_var(u, v) + " <= " + print_float(self.imp_bounds[(u, v)])
                self.edge_delay_bounds.append(bound)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def generate_td_imp_bounds(self, imp_bounds, tight):
        """Generates the improvement bounds.

        Parameters
        ----------
        imp_bounds : Dict[Tuple[str], Tuple[float]]
            Bounds on improvement for all edges. First element of the
            pair is the lower and the second the upper bound.
        tight : bool
            Decides if the bounds are enforced. If so, equality is assigned
            on the upper bound.
        
        Returns
        -------
        None
        """

        #Drop the previous bounds:
        self.edge_delay_bounds = []

        #String conversions:
        imp_var = lambda u, v : "imp_%d_%d" % (self.node_counts[u],\
                                               self.node_counts[v])
        rel = " = " if tight else " <= "

        for u, v in imp_bounds:
            bound = "" if tight else print_float(imp_bounds[(u, v)][0]) + " <= "\
                  + imp_var(u, v) + rel + print_float(imp_bounds[(u, v)][1])
            self.edge_delay_bounds.append(bound)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def generate_constraints(self):
       """Generates all constraints.
       Parameters
       ----------
       None
 
       Returns
       -------
       None
       """

       super(ImprovementGraph, self).generate_constraints()
       self.generate_td_imp_constraints()
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def print_constraints(self):
        """Prints all constraints.

        Parameters
        ----------
        None

        Returns
        -------
        str
            Text of the constraints.
        """

        txt = super(ImprovementGraph, self).print_constraints()

        txt += "\*   Edge Delays   *\\\n"
        for cst in self.edge_delay_constraints:
            txt += cst + "\n"

        return txt
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def print_bounds(self):
        """Prints all bounds.

        Parameters
        ----------
        None

        Returns
        -------
        str
            Text of the bounds.
        """

        txt = super(ImprovementGraph, self).print_bounds()
        for bound in self.edge_delay_bounds:
            if not bound:
                continue
            txt += bound + "\n"

        return txt
    #-----------------------------------------------------------------------#
##########################################################################

##########################################################################
class Selector(object):
    """Models the process of selecting the edges for improvement.

    Parameters
    ----------
    tg : networkx.DiGraph
        Timing graph with current edge delays annotated
        (using the >>td<< key). Cluster coordinates of each
        node are also annotated (using the >>coords<< key).
    timer : Timer
        A timer model, used for obtaining improvements.
    window : int
        The size of the window in which each node can move.

    Attributes
    ----------
    *Parameters
    relaxation : float
        The factor determining the amount of relaxation of the
        minimum achievable delay.
    fixed_nodes : List[str]
        Nodes that are to remain unmoved during ILP solving.
        For them, incident edge improvements amount to the shortest
        achievable programmable conenction, upon moving the other node.
    node_counts : Dict[str, int]
        Mapping between nodes and their counts.
    tg_lp : ImprovementGraph
        A timing graph model used for generation
        of arrival time and edge delay constraints.
    cov : Set[Tuple[str]]
        Set of edges that need to be improved.
    stretch : Set[Tuple[str]]
        Set of edges that are adjacent to the improved ones,
        and hence influenced by their improvement.
    max_delays : Dict[Tuple[str], float]
        Maximum allowed delays for all affected edges, such that
        the delay bound is still met.
        
    Methods
    -------
    set_relaxation(relaxation : float)
        Sets the relaxation factor.
    fix_nodes(List[str])
        Sets nodes that can not be moved by the ILP.
    find_corner_delay(u : str, v : str, min_max : str)
        Returns the minimum or the maximum delay of the given edge,
        given the current positions and fixing.
    find_bound()
        Returns the minimum achievable delay, under fixing constraints,
        but no improvement budget limitation.
    find_min_improvement()
        Finds some distribution of the minimum total positive
        improvement such that the relaxed bound is met.
    find_max_stretching()
        Finds some distribution of the maximum allowed delay increase
        of the edges adjacent to the improved ones, so that the
        relaxed bound is still met.
    """

    #-----------------------------------------------------------------------#
    def __init__(self, tg, timer, window):
        """Constructor of the Selector class.
        """
        #.......................................................................#
        def count_nodes(self):
            """Maps all nodes of the compressed graph to node counts.
    
            We assume that BLEs have already been merged and that
            FFs are designated by >>.d<< and >>.q<<
    
            Parameters
            ----------
            None
    
            Returns
            -------
            None
            """
   
            #FIXME: Duplicate from feeder.TimingGraph 
            self.node_counts = {}
            cnt = 0
            for u in sorted(tg.nodes()):
                if u.endswith(".d"):
                    rcnt = self.node_counts.get(u[:-1] + 'q', -1)
                    if rcnt >= 0:
                        self.node_counts.update({u : rcnt})
                        continue
                elif u.endswith(".q"):
                    rcnt = self.node_counts.get(u[:-1] + 'd', -1)
                    if rcnt >= 0:
                        self.node_counts.update({u : rcnt})
                        continue
                self.node_counts.update({u : cnt})
                cnt += 1
        #.......................................................................#

        self.tg = tg
        self.timer = timer
        self.window = window
        self.relaxation = 1
        self.fixed_nodes = []

        self.grid_width = max([c.x for c in self.timer.grid_types])
        self.grid_height = max([c.y for c in self.timer.grid_types])

        self.irrelevant = []

        #We want all edges to be improvable. If any is fixed by node
        #fixing, we can model that through zeroing the improvement bound.
        #Hence, we want all edges in ch and no nodes in tarrs,
        #apart from sources.
        count_nodes(self)
        self.inv_node_counts = {self.node_counts[u] : u[:-1] + 'd'\
                               if u.endswith(".q") else u\
                               for u in self.node_counts}

        E = [(u, v) for u, v in tg.edges()\
             if timer.grid_types[tg.node[u]["coords"]] == "clb"\
             and timer.grid_types[tg.node[v]["coords"]] == "clb"]
        self.ch = E

        tarrs = {u : 0 for u in tg if not tg.in_degree(u)}
        self.tg_lp = ImprovementGraph(self.tg, E, tarrs, self.node_counts) 
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def set_relaxation(self, relaxation):
        """Sets the relaxation factor.

        Parameters
        ----------
        relaxation : float
            The factor determining the amount of relaxation of the
            minimum achievable delay.

        Returns
        -------
        None

        Raises
        ------
        AssertionError
            If relaxation is not a float >= 1.
        """

        check_nonneg_float(relaxation)
        assert (relaxation >= 1.0), "Relaxation overconstraining."

        self.relaxation = relaxation
    #-----------------------------------------------------------------------#
    
    #-----------------------------------------------------------------------#
    def fix_nodes(nodes):
        """Sets nodes that can not be moved by the ILP.

        Parameters
        ----------
        nodes : List[str]
            List of nodes that should be fixed.

        Returns
        -------
        None

        Raises
        ------
        AssertionError
            If any of the nodes is not in the timing graph.
        """

        assert (all(self.tg.has_node(u) for u in nodes)),\
               "Trying to fix nodes that are not in the graph."

        self.fixed_nodes = nodes
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def find_corner_delay(self, u, v, min_max):
        """Returns the minimum or the maximum delay of the given edge,
        given the current positions and fixing.

        Parameters
        ----------
        u : str
            Tail node.
        v : str
            Head node.
        min_max : str
            Specifies whether the minimum or the maximum delay is being
            sought.

        Returns
        -------
        float
            The corner delay.
        """

        reverse = True if min_max == "max" else False

        #If both are fixed, no improvement is possible.
        if u in self.fixed_nodes and v in self.fixed_nodes:
            return self.tg[u][v]["td"]

        #If neither are fixed, we can pick the global minimum (maximum).
        #FIXME: In principle, we still want to bound endpoint movement,
        #which may mean that we are not going to strike the very global extreme.
        #NOTE: Fine for min, fix for max.
        if u not in self.fixed_nodes and v not in self.fixed_nodes:
            if min_max == "min":
                #FIXME: Almost duplicate from feeder.Stretcher.assign_cov
                min_td = self.tg[u][v]["td"]
                u_init_cls = self.tg.node[u]["coords"]
                v_init_cls = self.tg.node[v]["coords"]

                u_x_min = max(0, u_init_cls.x - self.window)
                u_x_max = min(self.grid_width, u_init_cls.x + self.window) + 1
                u_y_min = max(0, u_init_cls.y - self.window)
                u_y_max = min(self.grid_height, u_init_cls.y + self.window) + 1

                for x in range(u_x_min, u_x_max):
                    for y in range(u_y_min, u_y_max):
                        u_cls = NodeCls(x, y)
                        if self.timer.grid_types.get(u_cls, "<EMPTY>") != "clb":
                            continue
                        for dc in self.timer.direct_delays:
                            v_cls = NodeCls(u_cls.x + dc.x_offset, u_cls.y + dc.y_offset)
                            if self.timer.grid_types.get(v_cls, "<EMPTY>") == "clb":
                                cls_dc = PhysDirCon(NodePos(u_cls.x, u_cls.y, dc.u),\
                                                    NodePos(v_cls.x, v_cls.y, dc.v))
                                if abs(v_init_cls.x - v_cls.x) <= self.window\
                                    and abs(v_init_cls.y - v_cls.y) <= self.window:
                                    td = self.timer.time_cov((u, v), cls_dc)
                                    if td < min_td:
                                        min_td = td
                return min_td
            else:
                return max(self.timer.prog_delays["clb"]["clb"].values())\
                       * self.timer.scaleup


        u_cls = self.tg.node[u]["coords"]
        u_type = self.timer.grid_types[u_cls]
        v_cls = self.tg.node[v]["coords"]
        v_type = self.timer.grid_types[v_cls]

        #VPR differentiates only I/Os and others when constructing
        #the delay lookup tables.
        type_conv = lambda t : "io" if t == "io" else "clb"

        prog_delays = self.timer.prog_delays[type_conv(u_type)]\
                                            [type_conv(v_type)]
        prog_keys = sorted(prog_delays, key = lambda k : prog_delays[k],\
                           reverse = reverse)

        if u in self.fixed_nodes:
            for pk in prog_keys:
                new_v_cls_x = u_cls.x + pk.x_offset
                new_v_cls_y = u_cls.y + pk.y_offset
                if abs(v_cls.x - new_v_cls_x) <= self.window\
                   and abs(v_cls.y - new_v_cls_y) <= self.window:
                    con = PhysProgCon(u_cls, NodeCls(new_v_cls_x, new_v_cls_y))
                    td = self.timer.time_stretch((u, v), con)
                    if min_max == "min":
                        return min(td, self.tg[u][v]["td"])
                    else:
                        return max(td, self.tg[u][v]["td"])
        elif v in self.fixed_nodes:
            for pk in prog_keys:
                new_u_cls_x = v_cls.x - pk.x_offset
                new_u_cls_y = v_cls.y - pk.y_offset
                if abs(u_cls.x - new_u_cls_x) <= self.window\
                   and abs(u_cls.y - new_u_cls_y) <= self.window:
                    con = PhysProgCon(NodeCls(new_u_cls_x, new_u_cls_y), v_cls)
                    td = self.timer.time_stretch((u, v), con)
                    if min_max == "min":
                        return min(td, self.tg[u][v]["td"])
                    else:
                        return max(td, self.tg[u][v]["td"])
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def find_bound(self):
        """Returns the minimum achievable delay, under fixing
        constraints, but no improvement budget limitation.

        Parameters
        ----------
        None

        Returns
        -------
        float
            The minimum delay bound.
        """

        if not hasattr(self, "tg_lb"):
            self.tg_lb = self.tg.copy()

        for u, v in self.tg_lb.edges():
            self.tg_lb[u][v]["td"] = self.find_corner_delay(u, v, "min")

        cpd = nx.dag_longest_path_length(self.tg_lb, weight = "td")
        self.cpd = cpd

        return cpd
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def bound_improvements(self, E):
        """Returns a dictionary of improvement bounds for all edges in E.

        Parameters
        ----------
        None

        Returns
        -------
        Dict[Tuple[str], Tuple[float]]
            Dictionary of improvement bounds.
        """

        ff_strip = lambda u : u[:-len(".d")] if u.endswith((".d", ".q"))\
                              else u
        bound_dict = {}
        for u, v in E:
            if ff_strip(u) == ff_strip(v):
                ub = lb = self.tg[u][v]["td"]
            ub = self.tg[u][v]["td"] - self.tg_lb[u][v]["td"]
            lb = self.tg[u][v]["td"] - self.find_corner_delay(u, v, "max")
            bound_dict.update({(u, v) : (lb, ub)})
        
        return bound_dict
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def bound_pairwise_improvement(self, u):
        """For each pair of edges ((p, u), (u, c)) incident to the node u,
        finds the maximum combined improvement, given the allowed movement
        regions. FIXME: For now, only covering is considered.

        Parameters
        ----------
        u : str
            Node for which the pairs are considered.
        
        Returns
        -------
        List[str]
            A list of constraints.
        """

        #NOTE: Before compression, all LUT nodes are considered movable.
        #Hence, for any non-LUT neighbors, the improvement of the corresponding
        #edge will be fixed to 0, so we do not need to include such edges in
        #the pairwise computation here. Once we include stretching, we will
        #have to revisit this.

        ff_strip = lambda u : u[:-len(".d")] if u.endswith((".d", ".q"))\
                              else u

        parents = [p for p in self.tg.pred[u]\
                   if self.timer.grid_types[self.tg.node[p]["coords"]] == "clb"\
                   and ff_strip(p) != ff_strip(u) and (p, u) in self.ch]
        children = [c for c in self.tg[u]\
                   if self.timer.grid_types[self.tg.node[c]["coords"]] == "clb"\
                   and ff_strip(u) != ff_strip(c) and (u, c) in self.ch]

        if not hasattr(self, "in_offset_pairs"):
            #If offset pairs were not formed, do it now. This depends only on
            #the pattern, so it suffices to do it once.
            in_pairs = [(dc1, dc2, self.timer.direct_delays[dc1]\
                                 + self.timer.direct_delays[dc2])\
                                 for dc1 in self.timer.direct_delays\
                                 for dc2 in self.timer.direct_delays\
                        if dc1 != dc2 and dc1.v == dc2.v]
            in_pairs.sort(key = lambda p : p[2])
    
            in_pairs = [((-1 * p[0].x_offset, -1 * p[0].y_offset),\
                        (-1 * p[1].x_offset, -1 * p[1].y_offset),\
                        p[2]) for p in in_pairs]
     
            out_pairs = [(dc1, dc2, self.timer.direct_delays[dc1]\
                                 + self.timer.direct_delays[dc2])\
                                 for dc1 in self.timer.direct_delays\
                                 for dc2 in self.timer.direct_delays\
                        if dc1 != dc2 and dc1.u == dc2.u]
            out_pairs.sort(key = lambda p : p[2])
    
            out_pairs = [((p[0].x_offset, p[0].y_offset),\
                         (p[1].x_offset, p[1].y_offset),\
                         p[2]) for p in out_pairs]
          
            cross_pairs = [(dc1, dc2, self.timer.direct_delays[dc1]\
                                 + self.timer.direct_delays[dc2])\
                                 for dc1 in self.timer.direct_delays\
                                 for dc2 in self.timer.direct_delays\
                        if dc1 != dc2 and dc1.u == dc2.v]
            cross_pairs.sort(key = lambda p : p[2])
    
            cross_pairs = [((-1 * p[0].x_offset, -1 * p[0].y_offset),\
                           (p[1].x_offset, p[1].y_offset),\
                           p[2]) for p in cross_pairs]

            self.in_offset_pairs = in_pairs
            self.out_offset_pairs = out_pairs
            self.cross_offset_pairs = cross_pairs

        #String conversions:
        imp_var = lambda u, v : "imp_%d_%d" % (self.node_counts[u],\
                                               self.node_counts[v])
        bounds = []

        #.......................................................................#
        def pair_reachable(self, offset1, offset2, neighbor1, u, neighbor2):
            """Checks if there exists a movement of u and its two neighbors,
            so that they are aligned with the passed offsets.
            
            Parameters
            ----------
            offset1 : Tuple[int]
                Offset (as 2D vector) between u and the first neighbor.
            offset2 : Tuple[int]
                Offset (as 2D vector) between u and the second neighbor.
            neighbor1 : NodeCls
                Initial position of the first neighbor.
            u : NodeCls
                Initial position of u.
            neighbor2 : NodeCls
                Initial position of the second neighbor.

            Returns
            -------
            bool
                True if there is such a movement, False otherwise.
            """

            
            u_x_min = max(0, u.x - self.window)
            u_x_max = min(self.grid_width, u.x + self.window) + 1
            u_y_min = max(0, u.y - self.window)
            u_y_max = min(self.grid_height, u.y + self.window) + 1

            neighbor1_x_min = max(0, neighbor1.x - self.window)
            neighbor1_x_max = min(self.grid_width, neighbor1.x + self.window) + 1
            neighbor1_y_min = max(0, neighbor1.y - self.window)
            neighbor1_y_max = min(self.grid_height, neighbor1.y + self.window) + 1

            neighbor2_x_min = max(0, neighbor2.x - self.window)
            neighbor2_x_max = min(self.grid_width, neighbor2.x + self.window) + 1
            neighbor2_y_min = max(0, neighbor2.y - self.window)
            neighbor2_y_max = min(self.grid_height, neighbor2.y + self.window) + 1

            #FIXME: Can be made much more efficient by decoupling axes and
            #solving equations. Not critical for now.
            for x in range(u_x_min, u_x_max):
                neighbor1_x = u_x_min + offset1[0]
                if neighbor1_x < neighbor1_x_min or neighbor1_x > neighbor1_x_max:
                    continue
                neighbor2_x = u_x_min + offset2[0]
                if neighbor2_x < neighbor2_x_min or neighbor2_x > neighbor2_x_max:
                    continue
                for y in range(u_y_min, u_y_max):
                    neighbor1_y = u_y_min + offset1[0]
                    if neighbor1_y < neighbor1_y_min or neighbor1_y > neighbor1_y_max:
                        continue
                    neighbor2_y = u_y_min + offset2[0]
                    if neighbor2_y < neighbor2_y_min or neighbor2_y > neighbor2_y_max:
                        continue
                    return True

            return False
        #.......................................................................#

        for p1 in parents:
            for p2 in parents:
                if p1 == p2:
                    continue
                found = False
                td = self.tg[p1][u]["td"] + self.tg[p2][u]["td"]
                for op in self.in_offset_pairs:
                    if pair_reachable(self, op[0], op[1],\
                                      self.tg.node[p1]["coords"],
                                      self.tg.node[u]["coords"],
                                      self.tg.node[p2]["coords"]):
                        #FIXME: Use the timer routines!
                        new_td = self.timer.scaleup * (self.timer.tg[p1][u]["stripped_td"]\
                                                    + self.timer.tg[p2][u]["stripped_td"]\
                                                    + op[2])
                        bound = imp_var(p1, u) + " + " + imp_var(p2, u)\
                              + " <= " + print_float(td - new_td)
                        found = True
                        break
                if not found:
                    bound = imp_var(p1, u) + " + " + imp_var(p2, u)\
                          + " <= " + print_float(max(self.individual_imp_bounds[(p1, u)][1],\
                                                     self.individual_imp_bounds[(p2, u)][1]))
                bounds.append(bound)

        for c1 in children:
            for c2 in children:
                if c1 == c2:
                    continue
                found = False
                td = self.tg[u][c1]["td"] + self.tg[u][c2]["td"]
                for op in self.out_offset_pairs:
                    if pair_reachable(self, op[0], op[1],\
                                      self.tg.node[c1]["coords"],
                                      self.tg.node[u]["coords"],
                                      self.tg.node[c2]["coords"]):
                        #FIXME: Use the timer routines!
                        new_td = self.timer.scaleup * (self.timer.tg[u][c1]["stripped_td"]\
                                                    + self.timer.tg[u][c2]["stripped_td"]\
                                                    + op[2])
                        bound = imp_var(u, c1) + " + " + imp_var(u, c2)\
                              + " <= " + print_float(td - new_td)
                        found = True
                        break
                if not found:
                    bound = imp_var(u, c1) + " + " + imp_var(u, c2)\
                          + " <= " + print_float(max(self.individual_imp_bounds[(u, c1)][1],\
                                                     self.individual_imp_bounds[(u, c2)][1]))
                bounds.append(bound)

        for p in parents:
            for c in children:
                found = False
                td = self.tg[p][u]["td"] + self.tg[u][c]["td"]
                for op in self.cross_offset_pairs:
                    if pair_reachable(self, op[0], op[1],\
                                      self.tg.node[p]["coords"],
                                      self.tg.node[u]["coords"],
                                      self.tg.node[c]["coords"]):
                        #FIXME: Use the timer routines!
                        new_td = self.timer.scaleup * (self.timer.tg[p][u]["stripped_td"]\
                                                    + self.timer.tg[u][c]["stripped_td"]\
                                                    + op[2])
                        bound = imp_var(p, u) + " + " + imp_var(u, c)\
                              + " <= " + print_float(td - new_td)
                        bounds.append(bound)
                        found = True
                        break
                if not found:
                    bound = imp_var(p, u) + " + " + imp_var(u, c)\
                          + " <= " + print_float(max(self.individual_imp_bounds[(p, u)][1],\
                                                     self.individual_imp_bounds[(u, c)][1]))
                bounds.append(bound)

        return bounds
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def bound_degree_improvement(self, u):
        """Bounds the maximum in-, out-, and total-degree improvement of u.
        FIXME: For now this is a simplified version, not considering
        movement restrictions. This is not so bad, because the individual
        and pairwise bounds will somewhat constrict some of the addends.

        Parameters
        ----------
        u : str
            Node for which the degrees are considered.
        
        Returns
        -------
        List[str]
            A list of constraints.
        """

        indeg = self.tg.in_degree(u)
        outdeg = self.tg.out_degree(u)
        if indeg + outdeg <= 2:
            #These cases are treated by individal and pairwise bounds.
            return []

        #NOTE: Before compression, all LUT nodes are considered movable.
        #Hence, for any non-LUT neighbors, the improvement of the corresponding
        #edge will be fixed to 0, so we do not need to include such edges in
        #the computation here. Once we include stretching, we will have to revisit this.

        ff_strip = lambda u : u[:-len(".d")] if u.endswith((".d", ".q"))\
                              else u

        parents = [p for p in self.tg.pred[u]\
                   if self.timer.grid_types[self.tg.node[p]["coords"]] == "clb"\
                   and ff_strip(p) != ff_strip(u) and self.individual_imp_bounds[(p, u)][1]]
        children = [c for c in self.tg[u]\
                   if self.timer.grid_types[self.tg.node[c]["coords"]] == "clb"\
                   and ff_strip(u) != ff_strip(c) and self.individual_imp_bounds[(u, c)][1]]

        indeg = len(parents)
        outdeg = len(children)
        if indeg + outdeg <= 2:
            #Now strengthen the previous check.
            return []

        if indeg and not hasattr(self, "in_directs"):
            #If incoming directs were not formed, do it now. This depends only on
            #the pattern, so it suffices to do it once.
            in_directs = {}
            for dc in self.timer.direct_delays:
                try:
                    in_directs[dc.v].append((dc, self.timer.direct_delays[dc]))
                except:
                    in_directs.update({dc.v : [(dc, self.timer.direct_delays[dc])]})
            for pos in in_directs:
                in_directs[pos].sort(key = lambda p : p[1])
            self.in_directs = in_directs

        if outdeg and not hasattr(self, "out_directs"):
            #If outgoing directs were not formed, do it now. This depends only on
            #the pattern, so it suffices to do it once.
            out_directs = {}
            for dc in self.timer.direct_delays:
                try:
                    out_directs[dc.u].append((dc, self.timer.direct_delays[dc]))
                except:
                    out_directs.update({dc.u : [(dc, self.timer.direct_delays[dc])]})
            for pos in out_directs:
                out_directs[pos].sort(key = lambda p : p[1])
            self.out_directs = out_directs

        bounds = []

        in_directs = {}
        if indeg:
            in_directs = self.in_directs
        out_directs = {}
        if outdeg:
            out_directs = self.out_directs

        #String conversions:
        imp_var = lambda u, v : "imp_%d_%d" % (self.node_counts[u],\
                                               self.node_counts[v])

        #Time difference calculators:
        td_scaleup = self.timer.scaleup
        delta_in = lambda p : self.tg[p][u]["td"]\
                 - (self.timer.tg[p][u]["stripped_td"] * td_scaleup)
        delta_out = lambda c : self.tg[u][c]["td"]\
                 - (self.timer.tg[u][c]["stripped_td"] * td_scaleup)

        parent_deltas = sorted([(p, delta_in(p)) for p in parents], reverse = True)
        child_deltas = sorted([(c, delta_out(c)) for c in children], reverse = True)

        #.......................................................................#
        def direct_can_cover(self, u, v, dc):
            """Checks if the direct connection can be used to cover the edge (u, v(
    
            Parameters
            ----------
            u : str
                Tail node.
            v : str
                Head node.
            dc : ProgDirCon.
                Direct connection.
    
            Returns
            -------
            bool
                True if coverable, False otherwise.
            """
            
            u_init_cls = self.tg.node[u]["coords"]
            v_init_cls = self.tg.node[v]["coords"]
    
            u_x_min = max(0, u_init_cls.x - self.window)
            u_x_max = min(self.grid_width, u_init_cls.x + self.window) + 1
            u_y_min = max(0, u_init_cls.y - self.window)
            u_y_max = min(self.grid_height, u_init_cls.y + self.window) + 1
    
            for x in range(u_x_min, u_x_max):
                for y in range(u_y_min, u_y_max):
                    u_cls = NodeCls(x, y)
                    if self.timer.grid_types.get(u_cls, "<EMPTY>") != "clb":
                        continue
                    v_cls = NodeCls(u_cls.x + dc.x_offset, u_cls.y + dc.y_offset)
                    if self.timer.grid_types.get(v_cls, "<EMPTY>") == "clb":
                        if abs(v_init_cls.x - v_cls.x) <= self.window\
                           and abs(v_init_cls.y - v_cls.y) <= self.window:
                            return True
    
            return False
        #.......................................................................#

        max_in_imp = -1
        max_in_pos = None
        max_out_imp = -1
        max_out_pos = None
        max_total_imp = -1
        max_total_pos = None
        all_pos = list(in_directs.keys()) + list(out_directs.keys())
        for pos in all_pos:
            pos_in_directs = in_directs.get(pos, [])
            in_imp = 0
            if pos_in_directs:
                available_parents = copy.copy(parent_deltas)
                for dc in pos_in_directs:
                    matches = None
                    for p in available_parents:
                        if direct_can_cover(self, p[0], u, dc[0]):
                            matches = p
                            break
                    if matches is not None:
                        available_parents.remove(p)
                        in_imp += p[1] - (dc[1] * td_scaleup)
                        if not available_parents:
                            break
            if in_imp > max_in_imp:
                max_in_imp = in_imp
                max_in_pos = pos
            pos_out_directs = out_directs.get(pos, [])
            out_imp = 0
            if pos_out_directs:
                available_children = copy.copy(child_deltas)
                for dc in pos_out_directs:
                    matches = None
                    for c in available_children:
                        if direct_can_cover(self, u, c[0], dc[0]):
                            matches = c
                            break
                    if matches is not None:
                        available_children.remove(c)
                        out_imp += c[1] - (dc[1] * td_scaleup)
                        if not available_children:
                            break
            if out_imp > max_out_imp:
                max_out_imp = out_imp
                max_out_pos = pos
            total_imp = in_imp + out_imp
            if total_imp > max_total_imp:
                max_total_imp = total_imp
                max_total_pos = pos

        if max_in_pos is not None:
            bound = ""
            for p in parents:
                bound += imp_var(p, u) + " + "
            bound = bound[:-len(" + ")] + " <= " + print_float(max_in_imp)
            bounds.append(bound)
        if max_out_pos is not None:
            bound = ""
            for c in children:
                bound += imp_var(u, c) + " + "
            bound = bound[:-len(" + ")] + " <= " + print_float(max_out_imp)
            bounds.append(bound)
        if max_total_pos is not None:
            bound = ""
            for p in parents:
                bound += imp_var(p, u) + " + "
            for c in children:
                bound += imp_var(u, c) + " + "
            bound = bound[:-len(" + ")] + " <= " + print_float(max_total_imp)
            bounds.append(bound)

        return bounds
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def set_minimum_improvement(self, imp):
        """Sets the minimum improvement.
        This is useful to bound the minimum size of the core,
        thus maximizing the chances that it will still be successfully
        resolved.

        Parameters
        ----------
        imp : float
            Improvement bound.
        
        Returns
        -------
        None
        """

        self.min_total_imp = "total_imp >= " + print_float(imp) + "\n"
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def construct_improvement_lp(self):
        """Constructs and prints the improvement linear program.
        
        Parameters
        ----------
        None

        Returns
        -------
        str
            Text of the program in CPLEX format.
        """

        self.tg_lp.set_max_cpd(self.cpd * self.relaxation)
        self.tg_lp.generate_constraints()
        csts = self.tg_lp.print_constraints()

        bounds = self.bound_improvements(self.tg.edges())
        bounds = {b : (0, bounds[b][1]) for b in bounds}
        self.individual_imp_bounds = bounds

        nodes = [u for u in self.tg\
                 if self.timer.grid_types[self.tg.node[u]["coords"]] == "clb"]
        MAX_OUTDEG = 6
        for u in nodes:
            if self.tg.out_degree(u) <= MAX_OUTDEG:
                #NOTE: High-fanout nets can not be covered in any case.
                #For them, it is sufficient to bound the total improvement
                #over the whole incident-edge set. Pairwise would take too
                #much time even to generate.
                pairwise_csts = self.bound_pairwise_improvement(u)
                for pcst in pairwise_csts:
                    csts += pcst + "\n"
            degree_csts = self.bound_degree_improvement(u)
            for dcst in degree_csts:
                csts += dcst + "\n"
         
        self.tg_lp.generate_td_imp_bounds(bounds, tight = False)
        bounds = self.tg_lp.print_bounds()

        total_imp = "- total_imp + "
        for bound in self.tg_lp.edge_delay_bounds:
            edge_var = bound.split()[2]
            total_imp += edge_var + " + "
        csts += total_imp[:-len(" + ")] + " = 0\n"

        txt = "Minimize total_imp\n"
        txt += "Subject to\n"
        txt += csts

        txt += "Bounds\n"
        #We want to optionally set this externally.
        if hasattr(self, "min_total_imp"):
            bounds += self.min_total_imp
        txt += bounds
        txt += "End\n"

        return txt    
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def solve_lp(self, constructor, prob_name, timeout):
        """Calls the solver to solve the specified LP.

        Parameters
        ----------
        constructor : function
            The method that constructs the LP.
        prob_name : str
            Name of the problem, as used in filenames, without extension.
        timeout : int
            Time allowed for solving, in the number of seconds.

        Returns
        -------
        None

        Raises
        ------
        AssertionError
            If timeout is not a nonnegative integer.
        """

        check_nonneg_int(timeout)
        with open(prob_name + ".lp", "w") as prob:
            prob.write(constructor())

        if OPTIMIZER == "glpk":
            self.solve_glpk(prob_name, timeout)
        elif OPTIMIZER == "cplex":
            self.solve_cplex(prob_name, timeout)
        else:
            print("Unknown optimizer.")
            exit(-1)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def solve_cplex(self, prob_name, timeout):
        """Calls CPLEX to solve the LP.

        Parameters
        ----------
        prob_name : str
            Name of the problem, as used in filenames, without extension.
        timeout : int
            Time allowed for solving, in the number of seconds.

        Returns
        -------
        None
        """

        call = "cplex "
        call += "-c \"set timelimit %d\" " % timeout
        call += "\"read %s.lp\" \"optimize\" \"set logfile %s.sol\" "\
              % (prob_name, prob_name)
        call += "\"display solution variables -\" \"quit\""
        call += " > %s_cplex_call.log" % prob_name
        os.system(call)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def solve_glpk(self, prob_name, timeout):
        """Calls GLPK to solve the LP.

        Parameters
        ----------
        prob_name : str
            Name of the problem, as used in filenames, without extension.
        timeout : int
            Time allowed for solving, in the number of seconds.

        Returns
        -------
        None
        """

        call = "glpsol --tmlim %d --lp %s.lp" % (timeout, prob_name)
        call += " --write %s.sol --wglp %s.prob" % (prob_name, prob_name)
        call += " > %s.glpsol.log" % prob_name

        os.system(call)
    #-----------------------------------------------------------------------#
    
    #-----------------------------------------------------------------------#
    def parse_solution_cplex(self, prob_name):
        """Parses the CPLEX-produced solution, returning a dictionary of values
        of all variables.

        Parameters
        ----------
        None

        Returns
        -------
        Dict[str, float]
            Dictionary of variable values.
        """

        self.solution = {}
        try:
            with open("%s.sol" % prob_name, "r") as sol:
                lines = sol.readlines()
            #We must remove the solution file, or else it will be appended to
            #on the next run with the same problem name.
            os.system("rm -f %s.sol" % prob_name)
        except:
            return {}

        val_dict = {}
        rd = False
        for line in lines:
            if "Variable Name" in line:
                rd = True
                continue
            if not rd:
                continue
            if "All other variables" in line and "are 0" in line:
                break
            var = line.split()[0]
            if "**" in line.split()[1]:
                #This is how CPLEX designates values not within bounds,
                #in case the problem is infeasible.
                return {} 
            val = float(line.split()[1])
            val_dict.update({var : val})
        
        self.solution = val_dict

        return val_dict
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def parse_solution_glpk(self, prob_name):
        """Parses the GLPK-produced solution, returning a dictionary of values
        of all variables.

        Parameters
        ----------
        prob_name : str
            Name of the problem, as used in filenames, without extension.

        Returns
        -------
        Dict[str, float]
            Dictionary of variable values.
        """

        col_dict = {}
        with open(prob_name  + ".prob", "r") as prob:
            lines = prob.readlines()

        for line in lines:
            words = line.split()
            if len(words) != 4 or words[1] != 'j':
                continue
            col_dict.update({int(words[2]) : words[3]})

        with open(prob_name + ".sol", "r") as sol:
            lines = sol.readlines()

        val_dict = {}
        for line in lines:
            words = line.split()
            if len(words) != 5 or words[0] != 'j':
                continue
            var = col_dict[int(words[1])]
            val = float(words[3])
            val_dict.update({var : val})

        self.solution = val_dict

        return val_dict
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def parse_solution(self, prob_name):
        """Parses the solution, returning a dictionary of values of all
        variables.

        Parameters
        ----------
        prob_name : str
            Name of the problem, as used in filenames, without extension.

        Returns
        -------
        Dict[str, float]
            Dictionary of variable values.
        """

        if OPTIMIZER == "glpk":
            return self.parse_solution_glpk(prob_name)
        elif OPTIMIZER == "cplex":
            return self.parse_solution_cplex(prob_name)
        else:
            print("Unknown optimizer.")
            exit(-1)
    #-----------------------------------------------------------------------#
        
    #-----------------------------------------------------------------------#
    def get_improvements(self, min_imp):
        """Returns the improvements for all edges that should be improved.

        Parameters
        ----------
        min_imp : float
            Minimum improvement considered > 0.

        Returns
        -------
        Dict[Tuple[str], float]
            Minimum improvements required to achieve the relaxed bound.
        """
        
        #Fetchers
        fetch_u = lambda e : self.inv_node_counts[int(e.split("imp_")[1]\
                             .split('_')[0])]
        fetch_v = lambda e : self.inv_node_counts[int(e.split("imp_")[1]\
                             .split('_')[1])]

        minimum_improvements = {}
        for var in self.solution:
            try:
                u = fetch_u(var)
                if u.endswith(".d"):
                    u = u[:-1] + 'q'
                v = fetch_v(var)
                if v.endswith(".q"):
                    v = v[:-1] + 'd'
            except:
                continue
            val = self.solution[var]

            if abs(val) >= min_imp:
                minimum_improvements.update({(u, v) : val})

        return minimum_improvements
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def update_lb_graph_imp(self):
        """Updates the lower bound on delay of the selected edges
        with the achieved improvement after thresholding. Then
        runs STA again to produce the new critical path.

        Parameters
        ----------
        None

        Returns
        -------
        float
            The new critical path delay.
        """

        for u, v in self.tg_lb.edges():
            self.tg_lb[u][v]["td"] = self.tg[u][v]["td"]
        for e in self.min_edge_improvements:
            u, v = e
            self.tg_lb[u][v]["td"] -= self.min_edge_improvements[e]

        cpd = nx.dag_longest_path_length(self.tg_lb, weight = "td")
        self.cpd = cpd

        return cpd
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def select_edges_for_improvement(self, min_imp, timeout):
        """Selects the edges to be improved by constructing
        and solving the appropriate LP.

        Parameters
        ----------
        min_imp : float
            Minimum improvement considered > 0.
        timeout : int
            Time allowed for solving, in the number of seconds.

        Returns
        -------
        None
        """

        self.find_bound()
        self.solve_lp(self.construct_improvement_lp, "imp", timeout)
        self.parse_solution("imp")
        min_improvements = self.get_improvements(min_imp)

        self.min_edge_improvements = min_improvements
        affected_nodes = set([u for e in min_improvements for u in e])
        self.affected_nodes = affected_nodes.copy()
        for u in affected_nodes:
            if u.endswith(".d"):
                self.affected_nodes.add(u[:-1] + 'q')
            if u.endswith(".q"):
                self.affected_nodes.add(u[:-1] + 'd')

        self.update_lb_graph_imp()
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def construct_stretching_lp(self):
        """Constructs and prints the stretching linear program.
        
        Parameters
        ----------
        None

        Returns
        -------
        str
            Text of the program in CPLEX format.
        """

        #FIXME: Lose the hard code! (extra relaxation for thresholding)
        self.tg_lp.cpd += 10
        csts = self.tg_lp.print_constraints()

        #Fix the bounds of edges that are to be improved:
        bounds = {b : (0, self.min_edge_improvements[b])\
                      for b in self.min_edge_improvements}
        self.tg_lp.generate_td_imp_bounds(bounds, tight = True)
        #Find the stretchable edges:
        stretchable = set()
        for u in self.affected_nodes:
            for p in self.tg.pred[u]:
                if not (p, u) in self.min_edge_improvements:
                    stretchable.add((p, u))
            for c in self.tg[u]:
                if not (u, c) in self.min_edge_improvements:
                    stretchable.add((u, c))
        stretch_bounds = self.bound_improvements(stretchable)
        stretch_bounds = {s : (stretch_bounds[s][0], 0) for s in stretch_bounds}
        self.tg_lp.generate_td_imp_bounds(stretch_bounds, tight = False)

        ff_strip = lambda u : u[:-len(".d")] if u.endswith((".d", ".q"))\
                              else u
        self.stretchable = [k for k in stretch_bounds.keys()\
                            if ff_strip(k[0]) != ff_strip(k[1])]

        #Fix all other improvements to 0.
        other = {}
        for u, v in self.tg.edges():
            if not (u, v) in stretch_bounds\
               and not (u, v) in self.min_edge_improvements:
                other.update({(u, v) : (0, 0)})
        self.tg_lp.generate_td_imp_bounds(other, tight = True)

        bounds = self.tg_lp.print_bounds()

        txt = "Minimize "
        #We are still minimizing, because now the improvements are negative.

        #String conversions:
        imp_var = lambda u, v : "imp_%d_%d" % (self.node_counts[u],\
                                               self.node_counts[v])

        for u, v in stretch_bounds:
            txt += imp_var(u, v) + " + "
        txt = txt[:-len(" + ")] + "\n"

        txt += "Subject to\n"
        txt += csts

        txt += "Bounds\n"
        txt += bounds
        txt += "End\n"

        return txt    
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def maximize_stretching(self, timeout):
        """Finds a maximum total stretching of affected edges.

        Parameters
        ----------
        timeout : int
            Time allowed for solving, in the number of seconds.

        Returns
        -------
        None
        """
        self.solve_lp(self.construct_stretching_lp, "stretching", timeout)
        self.parse_solution("stretching")
        max_stretchings = self.get_improvements(-float("Inf"))
        #We want all stretchings, even if they are zero.

        max_delays = {}
        for e in self.stretchable:
            u, v = e
            max_delays.update({e : self.tg[u][v]["td"] - max_stretchings[e]})
            #Here improvements are negative, so we subtract them.

        self.max_delays = max_delays
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def prune_irrelevant_edges(self, max_achievable_delays, pruning_cpd = None):
        """Returns the list of connections that can not become critical
        under any possible movement of nodes. It does so by fixing the
        delays of all affected edges to their maximum under the movement
        region constraints (taken as input because the tightest regions
        are after conflict removal), then does an STA and return affected
        edges with nonnegative slack.

        Parameters
        ----------
        max_achievable_delays : Dict[Tuple[str], float]
            Maximum achievable delays under the given moving constraints
            and prior pruning.
        pruning_cpd : float
            Minimum expected CPD which is to be used for pruning irrelevant
            edges. If the target CPD is used for pruning, any objective
            value lower than the target may be invalid, because some edges
            that may invalidate it could have been removed.

        Returns
            List[Tuple[str]]
                A list of irrelevant edges.
        """

        if pruning_cpd is None:
            pruning_cpd = self.cpd

        tg_ub = self.tg.copy()
        for e in max_achievable_delays:
            u, v = e
            tg_ub[u][v]["td"] = max_achievable_delays[e]
        
        topo_nodes = list(nx.topological_sort(tg_ub))
        rev_topo_nodes = reversed(topo_nodes)
    
        for v in topo_nodes:
            tg_ub.node[v]["ta"] = max([tg_ub.node[u]["ta"] + tg_ub[u][v]["td"]\
                                    for u in tg_ub.pred[v]] + [0])
        for u in rev_topo_nodes:
            tg_ub.node[u]["tr"] = min([tg_ub.node[v]["tr"] - tg_ub[u][v]["td"]\
                                    for v in tg_ub[u]] + [pruning_cpd])

        total = 0
        irrelevant = []
        for u, v in tg_ub.edges():
            if not u in self.affected_nodes and not v in self.affected_nodes:
                continue
            slack = tg_ub.node[v]["tr"] - tg_ub.node[u]["ta"] - tg_ub[u][v]["td"]
            total += 1
            if slack >= 0:
                irrelevant.append((u, v))

        print "total", total

        self.irrelevant = irrelevant

        return irrelevant
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def sta_bound_delays(self):
        """Uses slack of each edge of the timing graph to determine its
        maximum tolerable delay. This way we still obtain some pruning, but
        we do not artificially create infeasible problems (e.g. p -3> u -7> c
        may be infeasible improvement, but p -5> u -5> c may be realistic; for
        the starting LP, they are the same of course, as it cares only of sums).
        When doing the STA we assume that each stretched edge keeps the initial
        delay, while each covered edge takes the lower bound delay.

        Parameters
        ----------
        None

        Returns
        -------
        Dict[Tuple[str], float]
            A dictionary of maximum delays for all edges.
        """

        tg_lb = self.tg.copy()
        for e in self.min_edge_improvements:
            u, v = e
            tg_lb[u][v]["td"] -= self.individual_imp_bounds[e][1]
 
        topo_nodes = list(nx.topological_sort(tg_lb))
        rev_topo_nodes = reversed(topo_nodes)
    
        for v in topo_nodes:
            tg_lb.node[v]["ta"] = max([tg_lb.node[u]["ta"] + tg_lb[u][v]["td"]\
                                    for u in tg_lb.pred[v]] + [0])
        for u in rev_topo_nodes:
            tg_lb.node[u]["tr"] = min([tg_lb.node[v]["tr"] - tg_lb[u][v]["td"]\
                                    for v in tg_lb[u]] + [self.cpd])
       
        max_tds = {}
        for u, v in tg_lb.edges():
            if not u in self.affected_nodes and not v in self.affected_nodes:
                continue
            max_td = tg_lb.node[v]["tr"] - tg_lb.node[u]["ta"]
            max_tds.update({(u, v) : max_td})

        return max_tds
    #-----------------------------------------------------------------------#
##########################################################################
