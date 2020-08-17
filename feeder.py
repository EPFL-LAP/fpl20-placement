"""Module that generates the inputs for ilp_placer
based on parsing of VPR files and results of improvement
edge selection process.

For the TimingGraph constructor, it takes a complete
timing graph as parsed from VPR's postplacement output,
along with improvement edge information. Runs the STA
and annotates the nodes accordingly. Then extracts the
convex subgraph on selected edges and adds constraining
dummy nodes. It assigns fixed arrival times to the dummy
parent nodes, all FF outputs, and the sink of the blocking
path (which models the longest path completely outside the
visible subgraph). All selected edges and all edges that
share a common node with them enter the changeable set.
Finally, it counts the BLEs to produce the node count map.
This is all that ilp_placer.TimingGraph needs.
"""

import networkx as nx
import copy
from config_vars import CROSSBAR_FEEDBACK_DELAY, CROSSBAR_INPUT_DELAY
from feeder_types import Ble, ble_empty
from ilp_placer_types import NodePos, NodeCls, BasicEdge,\
                             PhysDirCon, PhysProgCon, CovPair, StretchPair,\
                             FreeDirCon, FreeProgCon
from ilp_placer_types import node_pos_to_node_cls, phys_to_free_dir_con,\
                             phys_to_free_prog_con

NO_STRETCH_IMP = True

##########################################################################
class TimingGraph(object):
    """Class representing the timing graph and all the preprocessing
    manipulations that it requires.

    Parameters
    ----------
    tg : networkx.Digraph
        Timing graph annotated with edge delays (key = >>td<<).
    blif_nodes : List[str]
        Nodes of the timing graph that appear in the BLIF netlist.

    Attributes
    ----------
    graph : networkx.Digraph
        A copy of tg, used for manipulation.
    cov : List(Tuple(str))
        A list of edges selected for improvement by ILP.
    stretch : List(Tuple(str))
        A list of edges that were not selected for improvement
        but whose delay can change while ILP is being solved.
    movable : Set(str)
        A set of nodes movable by the ILP.
    fixed : Set(str)
        A set of nodes that are endpoint of edges whose delay can
        change, but are not allowed to be moved by the ILP.
    compressed_graph : networkx.Digraph
        A subgraph of graph, containing only the portion of the graph
        relevant for doing proper timing analysis (implicit) during
        ILP solving.
    node_counts : Dict[str, int]
        A mapping of nodes to natural numbers, with FF ports merged.

    Methods
    -------
    sta(tg : networkx.DiGraph)
        A vanilla STA routine.
    merge_bles(bles : List[Bles])
        Merges the BLEs together, by keeping only the FF.
    find_affected_edges(E : List[Tuple[str]])
        Given a set of edges for improving, finds all
        other edges affected by moving of their endpoint nodes.
    set_reachability(V : List[str])
        Performs reachability starting from the set V.
    set_is_reachable_from(u : str, V : List[str])
        Checks if any of the nodes of V is reachable from u.
    extract_convex_subgraph(E : List[Tuple[str]])
        Given a set of edges, extracts the convex subgraph.
    constrain_subgraph(subg : networkx.DiGraph)
        Adds constraining nodes to the subgraph with edge delays
        set to represent the missing portions of the original graph.
    add_blocking_path(subg : networkx.DiGraph)
        Adds a blocking path to a subgraph, modeling the longest
        path not in the subgraph.
    extact_contracted_subgraph(V : List[str])
        Given a set of nodes, extracts a convex subgraph on them,
        then contracts the paths between those nodes.
    preprocess_tg()
        Performs necessary preprocessing of the timing graph.
    compress_tg()
        Constructs the relevant subgraph of the timing graph.
    push_memory_dummies()
        Pushes the delay of the edges modeling the memory input
        set-up time to the edge driving the >>.d<< input and removes
        the appropriate dummy.
    count_nodes()
        Maps all nodes of the compressed graph to node counts.
    scale_delays(scaleup : float, tg : networkx.DiGraph)
        Scales all edge delays of the given graph.
    export_tg()
        Exports the compressed timing graph.
    export_tarrs()
        Exports arrival times of those nodes for which it can not
        change during ILP solving.
    export_node_counts()
        Exports node counts.  
    """

    #-----------------------------------------------------------------------#
    def __init__(self, tg, blif_nodes):
        """Constructor of the TimingGraph class.
        """

        self.graph = tg.copy()
        self.blif_nodes = blif_nodes
        self.preprocess_tg()
    #-----------------------------------------------------------------------#
    
    #-----------------------------------------------------------------------#
    def sta(self, tg):
        """A vanilla STA routine.
    
        Only annotates the arrival and required times of nodes.
        Required times of sinks are set to the length of the
        longest path, while the arrival times of sources are set
        to zero.
    
        Edge delay key is assumed to be "td" and arrival time
        is annotated with "ta" while required time is annotated
        with "tr".
    
        Parameters
        ----------
        tg : networkx.DiGraph
            Graph on which we want to run the STA.
 
        Returns
        -------
        Critical path delay.
        """

        cpd = nx.dag_longest_path_length(tg, weight = "td")
        topo_nodes = list(nx.topological_sort(tg))
        rev_topo_nodes = reversed(topo_nodes)
    
        for v in topo_nodes:
            tg.node[v]["ta"] = max([tg.node[u]["ta"] + tg[u][v]["td"]\
                                    for u in tg.pred[v]] + [0])
        for u in rev_topo_nodes:
            tg.node[u]["tr"] = min([tg.node[v]["tr"] - tg[u][v]["td"]\
                                    for v in tg[u]] + [cpd])
    
        return cpd
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def merge_bles(self, bles):
        """Merges the BLEs together, by keeping only the FF.

        Parameters
        ----------
        bles : List[Ble]
            A list of all BLEs in the packed netlist.

        Returns
        -------
        None
        """
        
        for ble in bles:
            if ble.lut and ble.ff:
                ff_d = ble.ff + ".d"
                ff_q = ble.ff + ".q"
                ble_td = self.graph[ble.lut][ff_d]["td"]
                #We assume that the Tsu and Tclk are already embedded
                #in the appropriate edge delays. LUT delay itself also.
                fanin = list(self.graph.pred[ble.lut])
                for u in fanin:
                    td = self.graph[u][ble.lut]["td"]
                    self.graph.add_edge(u, ff_d, td = td + ble_td)
                    self.graph.remove_edge(u, ble.lut)
                self.graph.remove_node(ble.lut)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def find_affected_edges(self, E):
        """Given the set of edges selected for improving, finds all
        other edges affected by moving of their endpoint nodes.

        E is stored as self.cov, while the remaining edges are stored
        as self.stretch.

        Movable and fixed affected nodes are stored as self.movable
        and self.fixed, respectively.

        Parameters
        ----------
        E : List[Tuple(str)]
            The set of edges selected for improvement.

        Returns
        -------
        None
        """

        self.cov = copy.deepcopy(E)
        self.stretch = []

        movable_nodes = set([u for e in E for u in e])
        missing_ff_pins = set()
        for u in movable_nodes:
            if u.endswith(".d"):
                missing_ff_pins.add(u[:-1] + 'q')
            if u.endswith(".q"):
                missing_ff_pins.add(u[:-1] + 'd')

        self.movable = movable_nodes | missing_ff_pins
        fixed_nodes = set([])

        ff_strip = lambda u : u[:-len(".d")] if u.endswith((".d", ".q"))\
                              else u

        for u in self.movable:
            for p in self.graph.pred[u]:
                if ff_strip(p) == ff_strip(u):
                    continue
                if not (p, u) in E:
                    self.stretch.append((p, u))
                if not p in self.movable:
                    fixed_nodes.add(p)
            for c in self.graph[u]:
                if ff_strip(u) == ff_strip(c):
                    continue
                if not (u, c) in E:
                    self.stretch.append((u, c))
                if not c in self.movable:
                    fixed_nodes.add(c)
        self.fixed = fixed_nodes
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def set_reachability(self, V):
        """Performs reachability starting from the set V.

        Parameters
        ----------
        V : List[str]
            A set of sources.
        
        Returns
        -------
        Set[str]
            All nodes reachable from V.
        """

        stack = list(V)
        reachable = set(V)
        while stack:
            v = stack.pop()
            for c in self.graph[v]:
                if not c in reachable:
                    stack.append(c)
                    reachable.add(c)
    
        return reachable
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def set_is_reachable_from(self, u, V):
        """Checks if any of the nodes of V is reachable from u.

        Parameters
        ----------
        u : str
            Node from which reachability is being checked.
        V : List[str]
            The set of target nodes.

        Returns
        -------
        bool
            True if a reachable node is found, False otherwise.
        """

        if u in V:
            return True

        stack = [u]
        visited = set([u])
        while stack:
            u = stack.pop()
            for v in self.graph[u]:
                if v in V:
                    return True
                if not v in visited:
                    stack.append(v)
                    visited.add(v)

        return False
    #-----------------------------------------------------------------------#
    
    #-----------------------------------------------------------------------#
    def extract_convex_subgraph(self, E):
        """Given a set of edges, extracts the convex subgraph.

        A set of nodes can also be passed.

        All nodes that are both reachable from heads of E and
        can reach tails of E are influenced by changes in E
        and should hence be part of the extracted subgraph.

        Parameters
        ----------
        E : List[Tuple[str]]
            Set of edges from which the subgraph is to be built.
            Alternatively, a flat list of string, representing
            nodes can be passed.

        Returns
        -------
        networkx.DiGraph
            The obtained convex subgraph.
        """

        if isinstance(E[0], str):
            affected = set(E)
        else:
            affected = set([u for e in E for u in e])

        reachable = self.set_reachability(affected)
        convex = affected.copy()
        for u in reachable:
            if self.set_is_reachable_from(u, affected):
                convex.add(u)
            
        return self.graph.subgraph(convex).copy()
    #-----------------------------------------------------------------------#
    
    #-----------------------------------------------------------------------#
    def constrain_subgraph(self, subg):
        """Adds constraining nodes to the subgraph with edge delays
        set to represent the missing portions of the original graph.

        Parameters
        ----------
        subg : networkx.DiGraph
            The subgraph to be constrained.
    
        Returns
        -------
        networkx.DiGraph
            Constrained subgraph.
        """

        sta_keys = set(["ta", "tr"])
        constrained = subg.copy()
        #.......................................................................#
        def next_dummy():
            """Finds the next unused dummy identifier.
            """
            dummy_cnt = 0
            while True:
                while constrained.has_node("dummy_%d" % dummy_cnt):
                    dummy_cnt += 1

                yield "dummy_%d" % dummy_cnt
        #.......................................................................#
        dummy = next_dummy()

        cpd = max([self.graph.node[u]["tr"] for u in self.graph])

        #Add periphery parents.
        for v in subg:
            for u in self.graph.pred[v]:
                if not subg.has_node(u):
                    d = next(dummy)
                    constrained.add_node(d, orig = u)
                    constrained.node[d].update({a : self.graph.node[u][a]\
                                                for a in self.graph.node[u]})
                    constrained.node[d]["ta"] = 0
                    constrained.node[d]["tr"] = cpd
                    td = self.graph.node[u]["ta"] + self.graph[u][v]["td"]
                    constrained.add_edge(d, v, td = td)
        #Add periphery children.
        for u in subg:
            for v in self.graph[u]:
                if not subg.has_node(v):
                    d = next(dummy)
                    constrained.add_node(d, orig = v)
                    constrained.node[d].update({a : self.graph.node[v][a]\
                                                for a in self.graph.node[v]})
                    constrained.node[d]["ta"] = 0
                    constrained.node[d]["tr"] = cpd
                    td = cpd - self.graph.node[v]["tr"] + self.graph[u][v]["td"]
                    constrained.add_edge(u, d, td = td)

        return constrained
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def add_blocking_path(self, subg):
        """Adds a blocking path to a subgraph, modeling the longest
        path not in the subgraph.

        Parameters
        ----------
        subg : networkx.DiGraph
            Subgraph to which the blocking path is to be added.

        Returns
        -------
        networkx.DiGraph
            Subgraph with the added blocking path.
        """

        cpd = nx.dag_longest_path_length(self.graph, weight = "td")
        blocked = subg.copy()
        complement = self.graph.subgraph([u for u in self.graph if not u in subg])
        external_cpd = nx.dag_longest_path_length(complement, weight = "td")

        dummy_cnt = 0
        while blocked.has_node("dummy_%d" % dummy_cnt):
            dummy_cnt += 1
        u = "dummy_%d" % dummy_cnt
        blocked.add_node(u, ta = 0, tr = cpd - external_cpd,\
                         coords = NodeCls(0, 0))
        while blocked.has_node("dummy_%d" % dummy_cnt):
            dummy_cnt += 1
        v = "dummy_%d" % dummy_cnt
        blocked.add_node(v, ta = external_cpd, tr = cpd,\
                         coords = NodeCls(0, 0))

        blocked.add_edge(u, v, td = external_cpd)

        return blocked
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def extract_contracted_subgraph(self, V):
        """Given a set of nodes, extracts a convex subgraph on them,
        then contracts the paths between those nodes.

        This is useful for stripping away the timing graph nodes
        repesenting many architectural features that have no corresponding
        BLIF nodes (e.g., BLE muxes).
        
        Parameters
        ----------
        V : List[str]
            The set of nodes on which the subgraph is constructed.

        Returns
        -------
        networkx.DiGraph
            The contracted convex subgraph.
        """

        self.sta(self.graph)
        subg = self.extract_convex_subgraph(V)
        nodes = list(subg.nodes())
        subg = self.constrain_subgraph(subg)
        subg = self.add_blocking_path(subg)
        constraint_nodes = [u for u in subg.nodes() if not u in nodes]

        relevant_nodes = set(V) | set(constraint_nodes)
        contracted = nx.DiGraph()
        contracted.add_nodes_from(relevant_nodes)
        for u in relevant_nodes:
            contracted.node[u].update(subg.node[u])

        visited = set([])
        for v in relevant_nodes:
            stack = [v]
            td_stack = [0]
            while stack:
                u = stack.pop(-1)
                visited.add(u)
                td = td_stack.pop(-1)
                for c in subg[u]:
                    edge_td = subg[u][c]["td"]
                    if c in relevant_nodes:
                        contracted.add_edge(v, c, td = td + edge_td)
                    else:
                        stack.append(c)
                        td_stack.append(td + edge_td)

        return contracted
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def preprocess_tg(self):
        """Performs the necessary preprocessing of the timing graph.

        After parsing, the timing graph contains a lot of nodes that
        model the architecture, but have no corresponding BLIF nodes.
        These need to be stripped. Likewise, all integer nodes are
        converted to strings, as expected here.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        #Convert node labels to strings.
        self.graph = nx.relabel_nodes(self.graph, {u : str(u)\
                                                   for u in self.graph})
        self.graph = self.extract_contracted_subgraph(self.blif_nodes)
    #-----------------------------------------------------------------------#
        
    #-----------------------------------------------------------------------#
    def compress_tg(self):
        """Constructs the relevant subgraph of the timing graph.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        ch = self.cov + self.stretch
        self.sta(self.graph)
        orig_edges = [(u, v) for u, v in self.graph.edges()]
        subg = self.extract_convex_subgraph(ch)
        subg = self.constrain_subgraph(subg)
        subg = self.add_blocking_path(subg)
        new_edges = [(u, v) for u, v in subg.edges if not (u, v) in orig_edges]

        for u, v in new_edges:
            if u in self.movable or v in self.movable:
                self.stretch.append((u, v))

        self.compressed_graph = subg
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def push_memory_dummies(self):
        """Pushes the delay of the edges modeling the memory input
        set-up time to the edge driving the >>.d<< input and removes
        the appropriate dummy. This way we can enforce the convention
        that no >>.d<< pin should have an outgoing edge.

        The changeable edge set is modified accordingly as well.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        edge_removal_list = []
        for u in self.compressed_graph:
            if u.startswith("memory") and u.endswith(".d"):
                c = list(self.compressed_graph[u])[0]
                p = list(self.compressed_graph.pred[u])[0]
                self.compressed_graph[p][u]["td"]\
                += self.compressed_graph[u][c]["td"]
                edge_removal_list.append((u, c))

        for u, v in edge_removal_list:
            self.compressed_graph.remove_edge(u, v)
            self.compressed_graph.remove_node(v)
            try:
                self.cov.remove((u, v))
            except:
                pass
            try:
                self.stretch.remove((u, v))
            except:
                pass
            try:
                self.movable.remove(v)
            except:
                pass
            try:
                self.fixed.remove(v)
            except:
                pass
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
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

        self.node_counts = {}
        cnt = 0
        for u in sorted(self.compressed_graph.nodes()):
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
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def scale_delays(self, scaleup, tg):
        """Scales all edge delays of the given graph.

        Parameters
        ----------
        scaleup : float
            Scaling factor (can be < 1.0)
        tg : networkx.DiGraph
            Graph on which the scaling is performed.        

        Returns
        -------
        None
        """
        
        for u, v in tg.edges():
            tg[u][v]["td"] *= scaleup
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def export_tg(self):
        """Exports the compressed timing graph.
        
        Parameters
        ----------
        None

        Returns
        -------
        networkx.DiGraph
            The compressed timing graph.
        """

        return self.compressed_graph
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def export_tarrs(self):
        """Exports arrival times of those nodes for which it can not
        change during ILP solving.

        This includes the FF outputs. We assume that FF ports are
        designated by >>.d<< and >>.q<<.

        Parameters
        ----------
        None

        Returns
        -------
        Dict[str, float]
            Constant arrival times.
        """

        #tg_lb = self.compressed_graph.copy()
        #for u, v in tg_lb.edges():
        #    if u in self.movable and v in self.movable:
        #        print u, v
        #        tg_lb[u][v]["td"] = 17.32
        #cpd = nx.dag_longest_path_length(tg_lb, weight = "td")
        #print cpd
        #print sorted(self.movable)
        #raw_input()

        tarrs = {}
        for u in nx.topological_sort(self.compressed_graph):
            if self.compressed_graph.in_degree(u) == 0:
                fix = True
            elif u in self.movable:
                fix = False
            else:
                parents = self.compressed_graph.pred[u]
                if all(p in tarrs for p in parents)\
                    and all(p not in self.movable for p in parents):
                    fix = True
                else:
                    fix = False
            if fix:
                tarrs.update({u : self.compressed_graph.node[u]["ta"]})

        return tarrs
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def export_node_counts(self):
        """Exports node counts.

        Parameters
        ----------
        None

        Returns
        -------
        Dict[str, int]
            Node counts.
        """

        return self.node_counts
    #-----------------------------------------------------------------------#
##########################################################################

##########################################################################
class Stretcher(object):
    """Models the process of deciding where each movable node can be
    placed, given a set of connections whose delay can change.

    Parameters
    ----------
    cov : Set[Tuple[str]]
        Set of edges whose delay should be improved.
    stretch : Set[Tuple[str]]
        Set of edges whose delay can change but need not improve.
    fixed_nodes : List[str]
        Set of nodes that should not be explicitly moved
        (they may be moved implicitly, but only in their own cluster).
    init_positions : Dict[str, NodeCls]
        Initial clusters of the movable nodes.
    dir_cons : List[FreeDirCon]
        A list of pattern connections.
    grid_types : Dict[NodeCls, str]
        Types of grid locations.
    window : int
        The size of the window in which each node can move.
    protected_positions : Dict[NodePos, str]
        A set of protected positions, each of which can only be
        used by the one indicated node.
    timer : Timer
        A covering and stretching delay mode.

    Attributes
    ----------
    *Parameters
    cls_dir_cons : Dict[NodeCls, List[PhysDirCon]]
        Direct conenctions grouped by the tail cluster.

    Methods
    -------
    in_window(u : str, pos : NodeCls)
        Checks if u is in the allowed window around
        its initial position
    cluster_dir_cons()
        Groups self.dir_cons based on the tail cluster and stores
        them in a dictionary indexed by the cluster coordinates.
    assign_cov(e : Tuple[str])
        Determines all direct connections that can be used
        to cover e given the current positions of its endpoints
        and allowed move windows. 
    assign_stretch(e : Tuple[str])
        Determines all cluster pairs that can be used
        to cover e given the current positions of its endpoints
        and allowed move windows.
    export_cov_map(node_counts : Dict[str, int])
        Exports the cov_map needed by ILP constructor, which
        is a mapping from edges of the timing graph whose delay
        can change during ILP solving and the physical direct
        connections of the FPGA grid.
    export_stretch_map(node_counts : Dict[str, int])
        Exports the stretch_map needed by ILP constructor, which
        is a mapping from edges of the timing graph whose delay
        can change during ILP solving and the pairs of clusters
        of the FPGA grid.
    """

    #-----------------------------------------------------------------------#
    def __init__(self, cov, stretch, fixed_nodes, init_positions,\
                 dir_cons, grid_types, window, protected_positions, timer):
        """Constructor of the Stretcher class.
        """

        self.cov = cov
        self.stretch = stretch
        self.fixed_nodes = fixed_nodes
        self.init_positions = init_positions
        self.dir_cons = dir_cons
        self.grid_types = grid_types
        self.window = {u : window for u in init_positions}
        self.protected_positions = protected_positions
        self.timer = timer

        self.valid_clusters = {}
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def override_window(self, u, window):
        """Overrides the default window of node >>u<< with >>window<<.

        Paramters
        ---------
        u : str
            Node for which the override should happen.
        window : int
            The new window size.
        
        Returns
        -------
        None
        """

        self.window[u] = window
    #-----------------------------------------------------------------------#
       
    #-----------------------------------------------------------------------#
    def in_window(self, u, pos):
        """Checks if u is in the allowed window around
           its initial position.
        
        Parameters
        ----------
        u : str
            Node for which window membership is being checked.
        pos : NodeCls
            New position for which window membership is being checked.
        
        Returns
        -------
        bool
            True if in window, False otherwise.
        """
        
        center = self.init_positions[u]

        if u in self.fixed_nodes:
            return pos == center

        if abs(center.x - pos.x) <= self.window[u]\
           and abs(center.y - pos.y) <= self.window[u]:
            return True

        return False
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def cluster_dir_cons(self):
        """Groups self.dir_cons based on the tail cluster and stores
        them in a dictionary indexed by the cluster coordinates.

        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """

        dir_con_dict = {}
        for u_cls in self.grid_types:
            if self.grid_types[u_cls] == "clb":
                for dc in self.dir_cons:
                    v_cls = NodeCls(u_cls.x + dc.x_offset, u_cls.y + dc.y_offset)
                    if self.grid_types.get(v_cls, "<EMPTY>") == "clb":
                        cls_dc = PhysDirCon(NodePos(u_cls.x, u_cls.y, dc.u),\
                                            NodePos(v_cls.x, v_cls.y, dc.v))
                        try:
                            dir_con_dict[u_cls].append(cls_dc)
                        except:
                            dir_con_dict.update({u_cls : [cls_dc]})

        self.cls_dir_cons = dir_con_dict
    #-----------------------------------------------------------------------#
 
    #-----------------------------------------------------------------------#
    def assign_cov(self, e):
        """Determines all direct connections that can be used
        to cover e given the current positions of its endpoints
        and allowed move windows. 

        Parameters
        ----------
        e : Tuple[str]
            The edge for which the covering pairs are being sought.

        Returns
        -------
        List[PhysDirCon]
            List of all valid direct connections that could cover e.
        """

        u, v = e
        covs = []

        u_init_cls = self.init_positions[u]
        v_init_cls = self.init_positions[v]
        u_type = self.grid_types[u_init_cls]
        v_type = self.grid_types[v_init_cls]
        if u_type != "clb" or v_type != "clb":
            return []

        if not u in self.valid_clusters:
            self.valid_clusters.update({u : set()})
        if not v in self.valid_clusters:
            self.valid_clusters.update({v : set()})

        grid_width = max([c.x for c in self.grid_types])
        grid_height = max([c.y for c in self.grid_types])

        u_x_min = max(0, u_init_cls.x - self.window[u])
        u_x_max = min(grid_width, u_init_cls.x + self.window[u]) + 1
        u_y_min = max(0, u_init_cls.y - self.window[u])
        u_y_max = min(grid_height, u_init_cls.y + self.window[u]) + 1

        if u in self.fixed_nodes:
            u_x_min = u_init_cls.x
            u_x_max = u_x_min + 1
            u_y_min = u_init_cls.y
            u_y_max = u_y_min + 1

        for x in range(u_x_min, u_x_max):
            for y in range(u_y_min, u_y_max):
                cls = NodeCls(x, y)
                if self.grid_types.get(cls, "<EMPTY>") != "clb":
                    continue
                for dc in self.cls_dir_cons.get(cls, []):
                    if self.protected_positions.get(dc.u, u) != u:
                        continue
                    if self.protected_positions.get(dc.v, v) != v:
                        continue 
                    if self.in_window(v, node_pos_to_node_cls(dc.v)):
                        td = self.timer.time_cov(e, dc)
                        covs.append(CovPair(dc.u, dc.v, td))
                        self.valid_clusters[u].add(node_pos_to_node_cls(dc.u))
                        self.valid_clusters[v].add(node_pos_to_node_cls(dc.v))

        return covs
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def assign_stretch(self, e):
        """Determines all cluster pairs that can be used
        to cover e given the current positions of its endpoints
        and allowed move windows.

        Parameters
        ----------
        e : Tuple[str]
            The edge for which the stretching pairs are being sought.

        Returns
        -------
        List[PhysProgCon]
            List of all valid cluster pairs that could cover e.
        """
        
        u, v = e
        stretches = []

        if not u in self.valid_clusters:
            self.valid_clusters.update({u : set()})
        if not v in self.valid_clusters:
            self.valid_clusters.update({v : set()})

        u_init_cls = self.init_positions[u]
        v_init_cls = self.init_positions[v]
        u_type = self.grid_types[u_init_cls]
        v_type = self.grid_types[v_init_cls]
        grid_width = max([c.x for c in self.grid_types])
        grid_height = max([c.y for c in self.grid_types])

        u_x_min = max(0, u_init_cls.x - self.window[u])
        u_x_max = min(grid_width, u_init_cls.x + self.window[u]) + 1
        u_y_min = max(0, u_init_cls.y - self.window[u])
        u_y_max = min(grid_height, u_init_cls.y + self.window[u]) + 1
        v_x_min = max(0, v_init_cls.x - self.window[u])
        v_x_max = min(grid_width, v_init_cls.x + self.window[u]) + 1
        v_y_min = max(0, v_init_cls.y - self.window[u])
        v_y_max = min(grid_height, v_init_cls.y + self.window[u]) + 1

        if u in self.fixed_nodes:
            u_x_min = u_init_cls.x
            u_x_max = u_x_min + 1
            u_y_min = u_init_cls.y
            u_y_max = u_y_min + 1
        if v in self.fixed_nodes:
            v_x_min = v_init_cls.x
            v_x_max = v_x_min + 1
            v_y_min = v_init_cls.y
            v_y_max = v_y_min + 1

        for u_x in range(u_x_min, u_x_max):
            for u_y in range(u_y_min, u_y_max):
                u_cls = NodeCls(u_x, u_y)
                if self.grid_types.get(u_cls, "<EMPTY>") != u_type:
                    continue
                for v_x in range(v_x_min, v_x_max):
                    for v_y in range(v_y_min, v_y_max):
                        v_cls = NodeCls(v_x, v_y)
                        if self.grid_types.get(v_cls, "<EMPTY>") != v_type:
                            continue
                        con = PhysProgCon(u_cls, v_cls)
                        td = self.timer.time_stretch(e, con)
                        if NO_STRETCH_IMP:
                            cur_con = PhysProgCon(u_init_cls, v_init_cls)
                            cur_td = self.timer.time_stretch(e, cur_con)
                            if cur_td > td:
                                td = cur_td
                        stretches.append(StretchPair(u_cls, v_cls, td))
                        self.valid_clusters[u].add(u_cls)
                        self.valid_clusters[v].add(v_cls)
       
        return stretches
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def export_cov_map(self, node_counts):
        """Exports the cov_map needed by ILP constructor, which
        is a mapping from edges of the timing graph whose delay
        can change during ILP solving and the physical direct
        connections of the FPGA grid.

        Parameters
        ----------
        node_counts : Dict[str, int]
            Mapping from node identifiers to nonnegative integers, used
            to mark variables in the ILP formulation.

        Returns
        -------
        cov_map : Dict[BasicEdge, List[CovPair]]
            A partial mapping of the timing-graph edge set to direct
            connections of the FPGA grid. Contains those edges whose
            delay is affected by the ILP and that can be covered by
            a direct connection.
        """

        cov_map = {}
        for e in self.cov:
            u, v = e
            covs = self.assign_cov(e)
            if covs:
                cov_map.update({BasicEdge(node_counts[u],\
                                          node_counts[v]) : covs})

        return cov_map
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def export_stretch_map(self, node_counts):
        """Exports the stretch_map needed by ILP constructor, which
        is a mapping from edges of the timing graph whose delay
        can change during ILP solving and the pairs of clusters
        of the FPGA grid.

        Parameters
        ----------
        node_counts : Dict[str, int]
            Mapping from node identifiers to nonnegative integers, used
            to mark variables in the ILP formulation.

        Returns
        -------
        stretch_map : Dict[BasicEdge, List[StretchPair]]
            A partial mapping of the timing-graph edge set to cluster
            pairs of the FPGA grid. Contains those edges whose delay is
            affected by the ILP.
        """

        stretch_map = {}
        for e in self.cov + self.stretch:
            stretches = self.assign_stretch(e)
            if stretches:
                u, v = e
                stretch_map.update({BasicEdge(node_counts[u],\
                                              node_counts[v]) : stretches})

        return stretch_map
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def resolve_node_conflicts(self, cov_map, stretch_map, inv_node_counts):
        """Constructs a valid position region for each node implied by
        each of the edges incident to it (depending on stretches and
        coverages) and prunes away all the pairs that imply a position
        outside the common region. If the common region reduces to zero,
        the problem has no feasible solution.
        """

        position_sets = {}
        cov_implied_positions = {}
        for e in set(cov_map.keys() + stretch_map.keys()):
            u, v = e
            covs = cov_map.get(e, [])
            stretches = stretch_map.get(e, [])
            u_clusters = [node_pos_to_node_cls(cov.u) for cov in covs]
            try:
                cov_implied_positions[u] |= set(u_clusters)
            except:
                cov_implied_positions.update({u : set(u_clusters)})
            v_clusters = [node_pos_to_node_cls(cov.v) for cov in covs]
            try:
                cov_implied_positions[v] |= set(v_clusters)
            except:
                cov_implied_positions.update({v : set(v_clusters)})
            u_clusters += [stretch.u for stretch in stretches]
            v_clusters += [stretch.v for stretch in stretches]

            u_clusters = set(u_clusters)
            v_clusters = set(v_clusters)
            try:
                position_sets[u].append(u_clusters)
            except:
                position_sets.update({u : [u_clusters]})
            try:
                position_sets[v].append(v_clusters)
            except:
                position_sets.update({v : [v_clusters]})


        for u in position_sets:
            if len(position_sets[u]) == 1:
                position_sets[u] = position_sets[u][0]
            else:
                position_sets[u] = set.intersection(*position_sets[u])
            try:
                u_init_cls = self.init_positions[inv_node_counts[u]]
            except:
                u_init_cls = self.init_positions[inv_node_counts[u][:-1] + 'q']
            if NO_STRETCH_IMP:
                position_sets[u] = [cls for cls in position_sets[u]\
                                    if cls in cov_implied_positions.get(u, [])\
                                    or cls == u_init_cls]
            if len(position_sets[u]) == 0:
                print(str(u) + ": Intersection empty. Problem infeasible")
            
       
        for e in cov_map:
            u, v = e
            rem_list = []
            for cov in cov_map[e]:
                u_cls = node_pos_to_node_cls(cov.u)
                if not u_cls in position_sets[u]:
                    rem_list.append(cov)
                else:
                    v_cls = node_pos_to_node_cls(cov.v)
                    if not v_cls in position_sets[v]:
                       rem_list.append(cov)
            for cov in rem_list:
                cov_map[e].remove(cov)
        for e in stretch_map:
            u, v = e
            rem_list = []
            for stretch in stretch_map[e]:
                u_cls = stretch.u
                if not u_cls in position_sets[u]:
                    rem_list.append(stretch)
                else:
                    v_cls = stretch.v
                    if not v_cls in position_sets[v]:
                       rem_list.append(stretch)
            for stretch in rem_list:
                stretch_map[e].remove(stretch)
                
        for u in position_sets:
            u_str = inv_node_counts[u]
            if not u_str in self.valid_clusters:
                if u_str.endswith(".d"):
                    u_str = u_str[:-1] + 'q'
                elif u_str.endswith(".q"):
                    u_str = u_str[:-1] + 'd'
            rem_list = []
            for cls in self.valid_clusters[u_str]:
                if not cls in position_sets[u]:
                    rem_list.append(cls)
            for cls in rem_list:
                self.valid_clusters[u_str].remove(cls)
 
        return cov_map, stretch_map
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def export_max_achievable_delays(self):
        """Exports the maximum achievable delays for affected conenctions.
        This is useful for removing irrelevant edges from the model. If we
        were to use stretch maps for that, some stretches could have been
        already pruned due to being too slow, resulting in incorrect
        annotation of the most-pessimistic timing graph.
        
        Parameters
        ----------
        None

        Returns
        -------
        Dict[Tuple[str], float]
            A dictionary of maximum achievable delays for all edges.
        """

        maads = {}
        for e in self.cov + self.stretch:
            u, v = e
            max_td = -1
            for u_cls in self.valid_clusters[u] | set([self.init_positions[u]]):
                for v_cls in self.valid_clusters[v] | set([self.init_positions[v]]):
                    con = PhysProgCon(u_cls, v_cls)
                    td = self.timer.time_stretch(e, con)
                    if max_td < td:
                        max_td = td
            maads.update({e : max_td})

        return maads
    #-----------------------------------------------------------------------#
##########################################################################

##########################################################################
class Timer(object):
    """Models timing of all connections, under different
    stretching and covering.

    Parameters
    ----------
    tg : networkx.DiGraph
        Postplacement timing graph cropped to the relevant BLIF nodes.
    direct_delays : Dict[FreeDirCon, float]
        Mapping between the direct connections and their delays.
    prog_delays : Dict[str, Dict[str, Dict[FreeProgCon, float]]]
        Mapping between the programmable connections and their delays.
        Connections are grouped by tail and head type, with types being
        >>clb<< and >>io<<, as in VTR-7.
    grid_types : Dict[NodeCls, str]
        Types of grid locations.

    Attributes
    ----------
    *Parameters
        
    Methods
    -------
    set_scaleup(scaleup : float)
        Sets the delay scaleup factor.
    lookup_delay(u_cls, v_cls):
        Fetches the delay from the VPR's table.
    strip_routing_delay()
        Strips all of the routing delay from eacg edge delay,
        apart from some necessary muxing (e.g., BLE output).
    time_cov(e : Tuple[str], cov : PhysDirCon)
        Returns the delay of the edge after the given covering.
    time_stretch(e : Tuple[str], stretch : PhysProgCon)
        Returns the delay of the edge after the given stretching.
    """

    #-----------------------------------------------------------------------#
    def __init__(self, tg, direct_delays, prog_delays, grid_types):
        """Constructor of the Timer class.
        """

        self.crossbar_feedback_delay = CROSSBAR_FEEDBACK_DELAY
        self.crossbar_input_delay = CROSSBAR_INPUT_DELAY

        self.tg = tg.copy()
        self.direct_delays = direct_delays
        self.prog_delays = prog_delays
        self.grid_types = grid_types
        self.scaleup = 1.0

        self.strip_routing_delay()
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def set_scaleup(self, scaleup):
        """Sets the delay scaleup factor.

        Parameters
        ----------
        scaleup : float
            Factor with which each returned delay will be multiplied.

        Returns
        -------
        None
        """

        self.scaleup = scaleup
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def lookup_delay(self, u_cls, v_cls):
        """Fetches the delay from the VPR's table.

        Parameters
        ----------
        u_cls : NodeCls
            Tail node cluster.
        v_cls : NodeCls
            Head node cluster.

        Returns
        -------
        float
            Delay.
        """

        #VPR differentiates only I/Os and others when constructing
        #the delay lookup tables.
        type_conv = lambda t : "io" if t == "io" else "clb"

        con = FreeProgCon(abs(v_cls.x - u_cls.x), abs(v_cls.y - u_cls.y))
        u_type = type_conv(self.grid_types[u_cls])
        v_type = type_conv(self.grid_types[v_cls])

        return self.prog_delays[u_type][v_type][con]
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def strip_routing_delay(self):
        """Strips all of the routing delay from eacg edge delay,
        apart from some necessary muxing (e.g., BLE output).

        For the results to be correct, the >>td<< labeled delay must
        correspond to the postplacement delay between the current
        coordinates of the edge's endpoints.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        for u, v, attrs in self.tg.edges(data = True):
            u_cls = self.tg.node[u]["coords"]
            v_cls = self.tg.node[v]["coords"]
            td = attrs["td"]
            if u_cls == v_cls:
                stripped_td = td - self.crossbar_feedback_delay
            else:
                stripped_td = td - self.crossbar_input_delay\
                            - self.lookup_delay(u_cls, v_cls)
            self.tg[u][v]["stripped_td"] = stripped_td
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def time_cov(self, e, cov):
        """Returns the delay of the edge after the given covering.

        Parameters
        ----------
        e : Tuple[str]
            Edge of the circuit graph for which delay is measured.
        cov : PhysDirCon
            Direct connection used for covering.
        """

        return self.scaleup * (self.tg[e[0]][e[1]]["stripped_td"]\
               + self.direct_delays[phys_to_free_dir_con(cov)])
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def time_stretch(self, e, stretch):
        """Returns the delay of the edge after the given stretching.

        Parameters
        ----------
        e : Tuple[str]
            Edge of the circuit graph for which delay is measured.
        stretch : PhysProgCon
            Cluster pair used for stretching.
        """

        if stretch.u == stretch.v:
            return self.scaleup * (self.tg[e[0]][e[1]]["stripped_td"]\
                   + self.crossbar_feedback_delay)
  
        return self.scaleup * (self.tg[e[0]][e[1]]["stripped_td"]\
               + self.lookup_delay(stretch.u, stretch.v)\
               + self.crossbar_input_delay)
    #-----------------------------------------------------------------------#
##########################################################################
