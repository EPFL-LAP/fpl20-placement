import networkx as nx
from collections import namedtuple
from feeder import Timer
from ilp_placer_types import *

from config_vars import N, I

CRIT_EXP = 8.0
T_UPDATE_FREQ = "none"
    
Path = namedtuple("Path", ["length", "cost", "nodes"])

DEBUG_CP = False
if DEBUG_CP:
    CP = []

##########################################################################
class Legalizer(object):
    """The class modeling the legalization process.

    Parameters
    ----------
    tg : nx.DiGraph
        The complete (i.e., not compressed, as we need all nodes)
        timing graph with all nodes annotated with their current
        cluster coordinates. The delays should correspond to the
        programmable connection delays implied by the current
        node positions.
    grid_types : Dict[NodeCls, str]
        Types of grid locations.
    coverage : Dict[Tuple[str], PhysDirCon]
        The coverage dictionary produced by the placer.
    direct_delays : Dict[FreeDirCon, float]
        Mapping between the direct connections and their delays.
    prog_delays : Dict[str, Dict[str, Dict[FreeProgCon, float]]]
        Mapping between the programmable connections and their delays.
        Connections are grouped by tail and head type, with types being
        >>clb<< and >>io<<, as in VTR-7.

    Attributes
    ----------
    *Parameters
    clusters : Dict[NodeCls, Cluster]
        A dictionary of cluster object, indexed by their coordinates.
    timer : feeder.Timer
        A timer object used to measure the influence of node movement.

    Methods
    -------
    refresh_timing()
        Updates the edge delays of tg, by considering the node
        positioning implied by coverage.
    compute_criticalities()
        Computes the exponentiated edge criticalities, used for
        movement weighing.
    build_clusters()
        Builds the cluster objects.
    construct_grid_graph()
        Constructs a grid graph with neighboring CLBs interconnected.
    """

    #-----------------------------------------------------------------------#
    def __init__(self, tg, grid_types, coverage, direct_delays, prog_delays):
        """Constructor of the Legalizer Class.
        """

        self.tg = tg
        self.grid_types = grid_types
        self.coverage = coverage
        self.direct_delays = direct_delays
        self.prog_delays = prog_delays

        self.clusters = {}
        self.timer = Timer(tg.copy(), direct_delays, prog_delays, grid_types)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def refresh_timing(self, verbose = False):
        """Updates the edge delays of tg, by considering the node
           positioning implied by coverage.

        Parameters
        ----------
        verbose : Optional[bool]
            Enables printing.
    
        Returns
        -------
        None
        """

        affected_nodes = {}
        for e in self.coverage:
            u, v = e
            u_cls = node_pos_to_node_cls(self.coverage[e].u)
            v_cls = node_pos_to_node_cls(self.coverage[e].v)
            affected_nodes.update({u : u_cls})
            affected_nodes.update({v : v_cls})
            if u.endswith(".q"):
                affected_nodes.update({u[:-1] + 'd' : u_cls})
            if v.endswith(".d"):
                affected_nodes.update({v[:-1] + 'q' : v_cls})
            self.tg[u][v]["td"] = self.timer.time_cov(e, self.coverage[e])

        visited = set(self.coverage.keys())
        for u in affected_nodes:
            u_cls = affected_nodes[u]
            for p in self.tg.pred[u]:
                e = (p, u)
                if e in visited:
                    continue
                p_cls = affected_nodes.get(p, None)
                #FIXME: Perhaps we should now consider unintentional edge covering.
                if p_cls is None:
                    p_cls = self.tg.node[p]["coords"]
                stretch = PhysProgCon(p_cls, u_cls)
                self.tg[p][u]["td"] = self.timer.time_stretch(e, stretch)
                visited.add(e)
            for c in self.tg[u]:
                e = (u, c)
                if e in visited:
                    continue
                c_cls = affected_nodes.get(c, None)
                #FIXME: Perhaps we should now consider unintentional edge covering.
                if c_cls is None:
                    c_cls = self.tg.node[c]["coords"]
                stretch = PhysProgCon(u_cls, c_cls)
                self.tg[u][c]["td"] = self.timer.time_stretch(e, stretch)
                visited.add(e)

        #Finally, we should update the coordinates:
        for u in affected_nodes:
            orig_cls = self.tg.node[u]["coords"]
            new_cls = affected_nodes[u]
            if verbose:
                if orig_cls != new_cls:
                    print(u + " : " + str(orig_cls) + " -> " + str(new_cls))
            self.tg.node[u]["coords"] = new_cls
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def compute_criticalities(self):
        """Computes the exponentiated edge criticalities,
           used for movement weighing.
    
        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        cpd = nx.dag_longest_path_length(self.tg, weight = "td")
        topo_nodes = list(nx.topological_sort(self.tg))
        rev_topo_nodes = reversed(topo_nodes)
    
        for v in topo_nodes:
            self.tg.node[v]["ta"] = max([self.tg.node[u]["ta"] + self.tg[u][v]["td"]\
                                        for u in self.tg.pred[v]] + [0])
        for u in rev_topo_nodes:
            self.tg.node[u]["tr"] = min([self.tg.node[v]["tr"] - self.tg[u][v]["td"]\
                                        for v in self.tg[u]] + [cpd])

        for u, v in self.tg.edges():
            slack = self.tg.node[v]["tr"] - self.tg.node[u]["ta"] - self.tg[u][v]["td"]
            crit = 1 - float(slack) / cpd
            self.tg[u][v]["t_weight"] = crit ** CRIT_EXP
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def build_clusters(self):
        """Builds the cluster objects.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
    
        #NOTE: We assume that the BLEs are merged in the timing graph.


        affected_nodes = {}
        for e in self.coverage:
            u, v  = e
            u_pos = self.coverage[e].u.i
            v_pos = self.coverage[e].v.i
            affected_nodes.update({u : u_pos})
            affected_nodes.update({v : v_pos})
            if u.endswith(".q"):
                affected_nodes.update({u[:-1] + 'd' : u_pos})
            if v.endswith(".d"):
                affected_nodes.update({v[:-1] + 'q' : v_pos})

        luts = {}
        for u in self.tg:
            cls = self.tg.node[u]["coords"]
            if self.grid_types[cls] != "clb":
                continue
            u_sink = u
            if u.endswith(".q"):
                u_sink = u[:-1] + 'd'
            ins = set(self.tg.pred[u_sink])
            u_src = u
            if u.endswith(".d"):
                u_src = u[:-1] + 'q'
            pos = affected_nodes.get(u, -1)
            lut = Lut(pos = pos, I = ins, o = u_src, fixed = pos != -1)
            try:
                luts[cls].append(lut)
            except:
                luts.update({cls : [lut]}) 
 
        for cls in luts:
            cluster = Cluster(cls)
            for lut in luts[cls]:
                cluster.add_lut(lut)
            self.clusters.update({cls : cluster})

        #Add the empty clusters as well.
        for loc in self.grid_types:
            if self.grid_types[loc] == "clb" and not loc in self.clusters:
                self.clusters.update({loc : Cluster(loc)})
    #-----------------------------------------------------------------------#
    
    #-----------------------------------------------------------------------#
    def construct_grid_graph(self):
        """Constructs a grid graph with neighboring CLBs interconnected.
        
        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.min_grid_x = 0
        self.min_grid_y = 0
        self.max_grid_x = max([cls.x for cls in self.grid_types])
        self.max_grid_y = max([cls.y for cls in self.grid_types])

        #String conversion:
        node_id = lambda cls : "%d_%d" % (cls.x, cls.y)

        grid_graph = nx.Graph()
        for x in range(self.min_grid_x, self.max_grid_x + 1):
            for y in range(self.min_grid_y, self.max_grid_y + 1):
                cls = NodeCls(x, y)
                if self.grid_types[cls] != "clb":
                    continue
                u = node_id(cls)
                for left in range(x - 1, self.min_grid_x - 1, -1):
                    left_cls = NodeCls(left, y)
                    if self.grid_types[left_cls] == "clb":
                        grid_graph.add_edge(u, node_id(left_cls), d = x - left)
                        break
                for right in range(x + 1, self.max_grid_x + 1):
                    right_cls = NodeCls(right, y)
                    if self.grid_types[right_cls] == "clb":
                        grid_graph.add_edge(u, node_id(right_cls), d = right - x)
                        break 
                for bottom in range(y - 1, self.min_grid_y - 1, -1):
                    bottom_cls = NodeCls(x, bottom)
                    if self.grid_types[bottom_cls] == "clb":
                        grid_graph.add_edge(u, node_id(bottom_cls), d = y - bottom)
                        break
                for up in range(y + 1, self.max_grid_y + 1):
                    up_cls = NodeCls(x, up)
                    if self.grid_types[up_cls] == "clb":
                        grid_graph.add_edge(u, node_id(up_cls), d = up - y)
                        break

        self.grid_graph = grid_graph
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def price_movement(self, lut, s_cls, t_cls):
        """Returns the cost of movement of LUT >>lut<< from the cluster
        >>s_cls<< to the cluster >>t_cls<< as a pair consisting of the magnitude
        of the worst negative new slack of the incident edges and the net change
        in sum of the criticality-weighted delays of incident edges, where the
        second is used for tie breaking.

        Parameters
        ----------
        lut : Lut
            The LUT to be priced.
        s_cls : NodeCls
            The original cluster.
        t_cls : NodeCls
            The new cluster.

        Returns
        -------
        (float, float)
            The cost of the movement.

        Raises
        ------
        ValueError
            If the LUT is fixed.
        """

        if lut.fixed:
            print("Trying to move a fixed LUT.")
            raise ValueError

        sink = lut.o
        if sink.endswith('.q'):
            sink = sink[:-1] + 'd'

        s_cost = 0
        t_cost = 0
        wns = 0
        for u in lut.I:
            e = (u, sink)
            t_weight = self.tg[u][sink]["t_weight"]
            s_td = self.tg[u][sink]["td"]
            s_cost += s_td * t_weight 
            u_cls = self.tg.node[u]["coords"]
            stretch = PhysProgCon(u_cls, t_cls)
            t_td = self.timer.time_stretch(e, stretch)
            t_cost += t_td * t_weight
            #Update the WNS.
            s_slack = self.tg[u][sink]["slack"]
            t_slack = s_slack - (t_td - s_td)
            if t_slack < wns:
                wns = t_slack

        u = lut.o
        for v in self.tg[u]:
            e = (u, v)
            t_weight = self.tg[u][v]["t_weight"]
            s_td = self.tg[u][v]["td"]
            s_cost += s_td * t_weight 
            v_cls = self.tg.node[v]["coords"]
            stretch = PhysProgCon(t_cls, v_cls)
            t_td = self.timer.time_stretch(e, stretch)
            t_cost += t_td * t_weight
            #Update the WNS.
            s_slack = self.tg[u][v]["slack"]
            t_slack = s_slack - (t_td - s_td)
            if t_slack < wns:
                wns = t_slack

        return (abs(wns), t_cost - s_cost)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def find_paths(self, src, d):
        """Finds the candidate paths for LUT movement, as in
        Darav et al. >>Multi-Commodity Flow-Based Spreading in Commercial
        Analytic Placer<<, FPGA'19

        Parameters
        ----------
        src : str
            The overflowed cluster we are listing paths for.
        d : int
            Maximum quadratic movement of a LUT.
    
        Returns
        -------
        List[Path]
            A list of paths of the grid graph nodes.

        Notes
        -----
        We change the approach slightly, in that the move pricing is based on
        the VPR-style net delay weighing; we allow any partially empty cluster
        to accept LUTs (they allow only fully empty ones);
        """

        #.......................................................................#
        def price_flow_edge(self, u, v):
            """Returns the cost of the least expensive move of a lut from
            the cluster represented in the grid graph by the node u to the
            cluster represented in the grid graph by the node v.

            Parameters
            ----------
            u : str
                source cluster node in the grid graph.
            v : str
                target cluster node in the grid graph.

            Returns
            -------
            float
                the cost of the cheapest move.
            """

            #string conversion:
            node_cls = lambda u : NodeCls(int(u.split('_')[0]), int(u.split('_')[1]))

            u_cls = node_cls(u)
            v_cls = node_cls(v)

            luts = [lut for lut in self.clusters[u_cls].luts if not lut.fixed]

            return min([self.price_movement(lut, u_cls, v_cls) for lut in luts])
        #.......................................................................#

        visited = set([src])
        partial_paths = [Path(length = 0, cost = [0, 0], nodes = [src])]
        complete_paths = []

        get_cls = lambda u : self.clusters[NodeCls(int(u.split('_')[0]),\
                                                   int(u.split('_')[1]))]
        get_fill = lambda u : len(get_cls(u).luts)
        get_supply = lambda u : get_cls(u).overflowed()
        get_demand = lambda u : N - get_fill(u) if N > get_fill(u) else 0
        
        supply = get_supply(src)
        total_demand = 0
        
        while partial_paths:
            path = partial_paths.pop()
            u = path.nodes[-1]
            for v in self.grid_graph[u]:
                if v in visited:
                    continue
                length = path.length + self.grid_graph[u][v]['d']
                if length > d:
                    continue
                cost = path.cost
                delta_cost = price_flow_edge(self, u, v)
                cost[0] += delta_cost[0]
                cost[1] += delta_cost[1]
                extended_path = Path(length, cost, path.nodes + [v])
                visited.add(v)
                demand = get_demand(v)
                if demand:
                    complete_paths.append(extended_path)
                    total_demand += demand
                if demand < N:
                    #If the cluster is completely empty, we can stop.
                    partial_paths.append(extended_path)
                if total_demand >= supply:
                    break
            if total_demand >= supply:
                break

        return sorted(complete_paths, key = lambda p : p.cost)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def lut_move(self, supply, paths):
        """Moves the LUTs along the found paths as in Darav et al, FPGA'19.

        Parameters
        ----------
        supply : int
            The amount of overflow in the source cluster.
        paths : List[Path]
            Paths assigned to the cluster being legalized.
        
        Returns
        -------
        None

        Notes
        -----
        We assume that the paths are already sorted in increasing cost.
        This function also covers parts of the main algorhtm of Darav et al.
        """

        #.......................................................................#
        def move_cheapest(self, u, v):
            """Moves the cheapest LUT from cluster u to cluster v.

            Parameters
            ----------
            u : str
                source cluster node in the grid graph.
            v : str
                target cluster node in the grid graph.

            Returns
            -------
            None
            """

            #string conversion:
            node_cls = lambda u : NodeCls(int(u.split('_')[0]), int(u.split('_')[1]))

            u_cls = node_cls(u)
            v_cls = node_cls(v)

            luts = [lut for lut in self.clusters[u_cls].luts if not lut.fixed]

            lut = min([lut for lut in luts],\
                      key = lambda l : self.price_movement(l, u_cls, v_cls))

            #print "moving", lut, self.price_movement(lut, u_cls, v_cls)

            self.clusters[u_cls].remove_lut(lut)
            self.clusters[v_cls].add_lut(lut)
        #.......................................................................#

        get_cls = lambda u : self.clusters[NodeCls(int(u.split('_')[0]),\
                                                   int(u.split('_')[1]))]
        get_fill = lambda u : len(get_cls(u).luts)
        get_demand = lambda u : N - get_fill(u) if N > get_fill(u) else 0

        for path in paths:
            if supply:
                v = path.nodes.pop()
                supply_dec = 0
                while path.nodes:
                    u = path.nodes.pop()
                    demand = get_demand(v)
                    #NOTE: It seems reasonable to saturate the demand
                    #before moving along the path as the paths were enumerated
                    #based on the sink demands. This can always be reverted
                    #when processing other overflowed clusters, although there
                    #will potentially be some timing impact.
                    for m in range(0, min(supply, demand)):
                        move_cheapest(self, u, v)
                        if T_UPDATE_FREQ == "move":
                            self.commit_to_tg()
                            print self.sta()
                            self.compute_criticalities()
                        supply_dec += 1
                    v = u
                supply -= supply_dec
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def legalize(self):
        """Performs cluster legalization according to Darav et al., FPGA'19.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.construct_grid_graph()

        ovfd = sorted([c for c in self.clusters if self.clusters[c].overflowed()],\
                       key = lambda cls : self.clusters[cls].overflowed())
        #NOTE: Most likely, we can make a set and then pop from it.
        i = -1
        cls_str = lambda cls : "%d_%d" % (cls.x, cls.y)
        print self.sta()
        while ovfd:
            i += 1
            d = (i + 1) ** 2
            for c in ovfd:
                paths = self.find_paths(cls_str(c), d)
                self.lut_move(self.clusters[c].overflowed(), paths)
            ovfd = sorted([c for c in self.clusters if self.clusters[c].overflowed()],\
                           key = lambda cls : self.clusters[cls].overflowed())
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def commit_to_tg(self):
        """Commits the changes in the clusters to the timing graph.
        
        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        coord_dict = {}
        for cls in self.clusters:
            for lut in self.clusters[cls].luts:
                coord_dict.update({lut.o : (cls, (lut.pos if lut.fixed else -1))})
                if lut.o.endswith('.q'):
                    coord_dict.update({lut.o[:-1] + 'd' : (cls,\
                                       (lut.pos if lut.fixed else -1))})

        for u, v in self.tg.edges():
            e = (u, v)
            u_cls = coord_dict.get(u, None)
            u_pos = -1
            if not u_cls is None:
                u_pos = u_cls[1]
                u_cls = u_cls[0]
            else:
                u_cls = self.tg.node[u]["coords"]
            v_cls = coord_dict.get(v, None)
            v_pos = -1
            if not v_cls is None:
                v_pos = v_cls[1]
                v_cls = v_cls[0]
            else:
                v_cls = self.tg.node[v]["coords"]
            if u_pos >= 0 and v_pos >= 0:
                cov = PhysDirCon(NodePos(u_cls.x, u_cls.y, u_pos),\
                                 NodePos(v_cls.x, v_cls.y, v_pos))
                dir_con = phys_to_free_dir_con(cov)
                if dir_con in self.timer.direct_delays:
                    self.tg[u][v]["td"] = self.timer.time_cov(e, cov)
                    continue
            stretch = PhysProgCon(u_cls, v_cls)
            self.tg[u][v]["td"] = self.timer.time_stretch(e, stretch)

        for u in coord_dict:
            self.tg.node[u]["coords"] = coord_dict[u][0]
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def sta(self):
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

        #FIXME: Duplicate from feeder.py.

        tg = self.tg

        cpd = nx.dag_longest_path_length(tg, weight = "td")
        
        if DEBUG_CP:
            changed = False
            global CP
            cp = nx.dag_longest_path(tg, weight = "td")
            new_cp = []
            for i in range(0, len(cp) - 1):
                u = cp[i]
                v = cp[i + 1]
                new_cp.append((u, v, tg[u][v]["td"]))
                if len(CP) > i:
                    if CP[i] != new_cp[-1]:
                        print CP[i], new_cp[-1]
                        changed = True
            new_cp.append((cp[-2], cp[-1], tg[cp[-2]][cp[-1]]["td"]))
            if len(CP) >= len(new_cp):
                if CP[len(new_cp) - 1] != new_cp[-1]:
                    print CP[len(new_cp) - 1], new_cp[-1]
                    changed = True
            CP = new_cp
            if changed:
                raw_input()

        topo_nodes = list(nx.topological_sort(tg))
        rev_topo_nodes = reversed(topo_nodes)
    
        for v in topo_nodes:
            tg.node[v]["ta"] = max([tg.node[u]["ta"] + tg[u][v]["td"]\
                                    for u in tg.pred[v]] + [0])
        for u in rev_topo_nodes:
            tg.node[u]["tr"] = min([tg.node[v]["tr"] - tg[u][v]["td"]\
                                    for v in tg[u]] + [cpd])

        for u, v in tg.edges():
            tg[u][v]["slack"] = tg.node[v]["tr"] - tg.node[u]["ta"] - tg[u][v]["td"]

        return cpd
    #-----------------------------------------------------------------------#
##########################################################################

Lut = namedtuple("Lut", ["pos", 'I', 'o', "fixed"])
"""A named tuple describing a LUT.
   pos : int
       Position of the LUT within the cluster.
       -1 means that the exact position is not assigned.
   I : Set[str]
       The set of inputs of the LUT.
   o : str
       The output that the LUT provides.
   fixed : bool
       Tells if the LUT can be moved or not.
"""

##########################################################################
class Cluster(object):
    """Class modeling a Cluster object.

    Parameters
    ----------
    coords : NodeCls
        Coordinates of the cluster.

    Attributes
    ----------
    luts : List[Lut]
        LUTs that the cluster holds.
    inputs : Set[str]
        A set of inputs to the cluster.
    outputs : Set[str]
        A set of outputs of the cluster.

    Methods
    -------
    add_lut(lut : Lut)
        Adds a LUT to the cluster. Input and output lists are maintained.
    remove_lut(lut : LUT)
        Removes a LUT from the cluster. Input and output lists are maintained.
    overflowed()
        Returns the number of LUTs beyond the capacity.
    """

    #-----------------------------------------------------------------------#
    def __init__(self, coords):
        """Constructor of the Cluster class.
        """

        self.coords = coords
        self.luts = []
        self.inputs = set()
        self.outputs = set()
    #-----------------------------------------------------------------------#
    
    #-----------------------------------------------------------------------#
    def add_lut(self, lut):
        """Adds a LUT to the cluster. Input and output lists are maintained.
        
        Parameters
        ----------
        lut : Lut
            The lut to be added.

        Returns
        -------
        None
        """

        #FIXME: For now direct connections providing inputs are not considered.

        if lut in self.luts:
            return

        self.luts.append(lut)
        self.outputs.add(lut.o)
        self.inputs |= lut.I - self.outputs
        try:
            self.inputs.remove(lut.o)
        except:
            pass
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def remove_lut(self, lut):
        """Removes a LUT from the cluster. Input and output lists are maintained.
        
        Parameters
        ----------
        lut : Lut
            The lut to be added.

        Returns
        -------
        None

        Notes
        -----
        Attempting to remove a fixed LUT will be accepted silently, as this
        can still be useful for restructuring the cluster. Make sure that
        you are aware of such attempts.
        """

        #FIXME: For now direct connections providing inputs are not considered.

        try:
            self.luts.remove(lut)
            self.outputs.remove(lut.o)
        except:
            return

        #Check if the output is now required.
        o_req = False
        for cmp_lut in self.luts:
            if lut.o in cmp_lut.I:
                o_req = True
                break
        if o_req:
            self.inputs.add(lut.o)
        
        #Check if some inputs are no longer required.
        rem_list = []
        for i in lut.I:
            needed = False
            for cmp_lut in self.luts:
                if i in cmp_lut.I or i == cmp_lut.o:
                    needed = True
                    break
            if not needed:
                rem_list.append(i)
        #for i in rem_list:
        #    self.inputs.remove(i)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def overflowed(self):
        """Returns the number of LUTs beyond the capacity.

        Paramters
        ---------
        None

        Returns
        -------
        int
            The amount of overflow.
        """

        ovf = len(self.luts) - N
        if ovf < 0:
            ovf = 0

        return ovf
    #-----------------------------------------------------------------------#
##########################################################################
