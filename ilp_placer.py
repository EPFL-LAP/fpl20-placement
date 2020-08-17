"""This module contains everything that is necessary to set up an ILP
problem for placement in the presence of direct connections, solve it,
and extract the placement from the solution.
"""

import networkx as nx
import copy
import os
from ilp_placer_types import *
from ilp_placer_base import *

from config_vars import N, I

OPTIMIZER = "cplex"

FILTER_SLOW = True
#Remove cover and stretch pairs that are too slow.

ENFORCE_LEGALIZATION = False
#Enforce cluster size and input count legality.

DELAY_PIVOT = False
#Merge various pairs based on their delays.

LAZY_STRETCHES = False
#Move the stretching constraints to lazy.

TIGHT_TD_BOUNDS = True
#Bound edge welays based on actually achievable delays.

TIGHT_TA_BOUNDS = True
#Bound arrival times based on optimistic and pessimistic STA.

REMOVE_IRRELEVANT_EDGES = True
#Remove edges from the timing graph that can never become critical.

INCIDENCE_COUNTING = True
#Add fanin/fanout counting constraints and implications
#on node positions they create.

ONLY_COV_CLS = True
#If no covering pair dictates the cluster of a node, force it to the initial one.

MINIMIZE_DISPLACEMENT = False
#If feasibility is sought, sets the objective to minimum total quadratic cluster displacement.

MAXIMIZE_COVERAGE_SIMILARITY = False
#If feasibility is sought, sets the maximum similarity between the covered indicator
#set of the incumbent solution (passed externally) as the objective.

FORCE_COVERAGE_INCREASE = False
#Force each solution to have wider coverage (in terms of cardinality) than the incumbent.

OPT_MODE = "feas"
precision_4 = lambda f : "%.4f" % f
precision_0 = lambda f : "%d" % f
print_float = precision_0

##########################################################################
class Edge(object):
    """Class describing an edge of the timing graph which should
    enter the ILP formulation.

    Parameters
    ----------
    u : int
        Tail node.
    v : int
        Head node.
    init_td : float
        Initial delay (i.e. postplacement)

    Attributes
    ----------
    *Parameters    
    td_max : float
        Maximum tolerable delay.
    cov_pairs : List[CovPair]
        A list of position pairs that cover the given edge with a
        direct connection of a certain delay.
    stretch_pairs : List[StretchPair]
        A list of position pairs, in terms of clusters, that stretch
        the given conenction without covering it. Note that
        ``stratching'' can also reduce to translation or contraction. 
        
    Methods
    -------
    add_cov(u_pos : NodePos, v_pos : NodePos, td : float)
        Adds a covering position pair.
    add_stretch(u_pos : NodeCls, v_pos : NodeCls, td : float)
        Adds a noncovering position pair.
    set_max_td(td : float)
        Sets the maximum tolerable delay for the edge.
    prune()
        Prunes away all the endpoint position pairs that would result
        in exceeding the timing budget.
    set_tail_alias(alias : str)
        Sets the alias of the tail node.
    set_head_alias(alias : str)
        Sets the alias of the head node.
    generate_constraints(fixed_nodes : Dict[int, NodeCls])
        Generates constraints for the edge.
    print_constraints()
        Prints all constraints for the given edge in
        a single block.
    print_td_bound()
        Prints the delay bound.
    """

    #-----------------------------------------------------------------------#
    def __init__(self, u, v, init_td):
        """Constructor of the Edge class.
        """

        e = BasicEdge(u = u, v = v)
 
        self.u = u
        self.u_alias = str(u)
        self.v = v
        self.v_alias = str(v)
        self.td_max = float("inf")
        self.init_td = init_td
        self.cov_pairs = []
        self.stretch_pairs = []
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def add_cov(self, u_pos, v_pos, td):
        """Adds a covering position pair.

        Parameters
        ----------
        u_pos : NodePos
            Position of the tail node.
        v_pos : NodePos
            Position of the head node.
        td : float
            Delay of the edge if covered by the direct connection
            between u_pos and v_pos.
        
        Returns
        -------
        None
        """

        cp = CovPair(u = u_pos, v = v_pos, td = td)
        if FILTER_SLOW and td > self.td_max:
            return
        self.cov_pairs.append(cp)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def add_stretch(self, u_pos, v_pos, td):
        """Adds a noncovering position pair.

        Parameters
        ----------
        u_pos : NodeCls
            Position of the tail node.
        v_pos : NodeCls
            Position of the head node.
        td : float
            Delay of the edge if u and v are at u_pos and v_pos
            cluster, respectively.
        
        Returns
        -------
        None
        """

        sp = StretchPair(u = u_pos, v = v_pos, td = td)
        if FILTER_SLOW and td > self.td_max:
            return
        self.stretch_pairs.append(sp)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def set_max_td(self, td):
        """Sets the maximum tolerable delay for the edge.

        Parameters
        ----------
        td : float
            The maximum tolerable delay.
        
        Returns
        -------
        None
        """

        self.td_max = td
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def set_tail_alias(self, alias):
        """Sets the alias of the tail node.

            Useful for printing the ILP in a way that eases
            manual debug.

        Parameters
        ----------
        alias : str
            Alias of the tail node.
        
        Returns
        -------
        None
        """

        self.u_alias = alias
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def set_head_alias(self, alias):
        """Sets the alias of the head node.

            Useful for printing the ILP in a way that eases
            manual debug.

        Parameters
        ----------
        alias : str
            Alias of the head node.
        
        Returns
        -------
        None
        """

        self.v_alias = alias
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def prune(self):
        """Prunes away all the endpoint position pairs that would result
        in exceeding the timing budget.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        rem_list = []
        for i, cp in enumerate(self.cov_pairs):
            if FILTER_SLOW and cp.td > self.td_max:
                rem_list.append(i)
        for i in reversed(rem_list):
            self.cov_pairs.pop(i)

        rem_list = []
        for i, sp in enumerate(self.stretch_pairs):
            if FILTER_SLOW and sp.td > self.td_max:
                rem_list.append(i)
        for i in reversed(rem_list):
            self.stretch_pairs.pop(i)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def group_delays(self, fixed_nodes):
        """Generates the delay assignment constraint of the form:
        td = Sum(Tau * indicator), for all achievable constant delays Tau.
        Then indicator = Sum(pair), for all pairs that achieve the indicated
        delay. This way, the solver should be able to quickly eliminate large
        portions of the search space, by profiting from the timing constraints.

        Paramters
        ---------
        fixed_nodes : Dict[int, NodeCls]
        A mapping of nodes not movable by the ILP to their
        respective cluster positions. 

        Returns
        -------
        None
        """

        #String conversions:

        u = str(self.u)
        v = str(self.v)
        e = "_%d_%d" % (self.u, self.v)
        
        pos_suffix = lambda pos : "_%d_%d_%d" % (pos.x, pos.y, pos.i)
        cls_suffix = lambda cls : "_%d_%d" % (cls.x, cls.y)

        u_pos_var = lambda cp : "u_" + u + pos_suffix(cp.u)
        v_pos_var = lambda cp : "u_" + v + pos_suffix(cp.v)
        cp_var = lambda cp : 'e' + e + pos_suffix(cp.u) + pos_suffix(cp.v)

        u_cls_var = lambda stc : "u_" + u + cls_suffix(stc.u)
        v_cls_var = lambda stc : "u_" + v + cls_suffix(stc.v)
        stc_var = lambda stc : 'e' + e + cls_suffix(stc.u) + cls_suffix(stc.v)

        tau_var = lambda tau : "tau" + e + '_' + tau

        #Constraint declaration:

        td_assign = "- td" + e + " + "

        cov_taus = set()
        #NOTE: By doing some manipulation of delays (e.g., adding a picosecond),
        #      we can make all cov taus differ from any stretch tau. This way,
        #      we can compute the covering indices based on taus directly, thus
        #      making propagation potentially faster. Also, branching order can
        #      be changed (we do not need to know exactly how an edge is covered
        #      before we know if edges selected for covering can be covered at all).
        delay_classes = {}
        for cp in self.cov_pairs:
            td = print_float(cp.td)
            #NOTE: We want only the precision that will be seen by the solver.
            try:
                delay_classes[td].append((cp, "cov"))
            except:
                delay_classes.update({td : [(cp, "cov")]})
            cov_taus.add(td)
        for stc in self.stretch_pairs:
            td = print_float(stc.td)
            try:
                delay_classes[td].append((stc, "stretch"))
            except:
                delay_classes.update({td : [(stc, "stretch")]})

        if any([p for p in delay_classes[tau] if p[1] == "cov"]\
           and [p for p in delay_classes[tau] if p[1] == "stretch"]\
           for tau in delay_classes):
            print "Coverage not derivable from Taus."
            exit(-1)
 
        self.td_assign = ""
        if not hasattr(self, "bin_vars"):
            self.bin_vars = set()
        lin_prefix = 'u' if self.cov_pairs else ''

        tau_force = ""
        for tau in delay_classes:
            tau_assignment = "- " + tau_var(tau) + " + "
            tau_force += tau_var(tau) + " + "
            self.bin_vars.add(tau_var(tau))
            for pair in delay_classes[tau]:
                if pair[1] == "cov":
                    cp = pair[0]
                    extend_var = cp_var(cp)
                    if self.u in fixed_nodes:
                        extend_var = v_pos_var(cp)
                    elif self.v in fixed_nodes:
                        extend_var = u_pos_var(cp)
                else:
                    stc = pair[0]
                    extend_var = lin_prefix + stc_var(stc)
            
                #Extend tau assingment:
                tau_assignment += extend_var + " + "
            tau_assignment = tau_assignment[:-len(" + ")] + " = 0\n"
            self.td_assign += tau_assignment

            #Extend the delay assignment constraint:
            td_assign += print_float(pair[0].td) + ' ' + tau_var(tau) + " + "
        td_assign = td_assign[:-len(" + ")] + " = 0\n"
        self.td_assign += td_assign
        tau_force = tau_force[:-len(" + ")] + " = 1\n"
        self.td_assign += tau_force

        self.cov_assign = ""
        if self.cov_pairs:
            for tau in cov_taus:
                self.cov_assign += tau_var(tau) + " + "
            self.cov_assign = self.cov_assign[:-len(" + ")]\
                            + " - cov_ind" + e + " = 0"
            self.bin_vars.add("cov_ind" + e)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def generate_constraints(self, fixed_nodes):
        """Generates constraints for the edge.

        The following types of constraints are generated:

        1) Linearizations of all endpoint position products.
        2) Covering indicator constraint.
        3) Uncovering indicator constraint, if there are both
           covering and stretching pairs.
        4) Delay assignment.
        5) Linearizations of stretching-pair delays, if conditions
           of 3) are satisfied.
        6) Delay bound, depending on the given td_max.

        Paramters
        ---------
        fixed_nodes : Dict[int, NodeCls]
        A mapping of nodes not movable by the ILP to their
        respective cluster positions. 

        Returns
        -------
        None
        """

        #String conversions:

        u = str(self.u)
        v = str(self.v)
        e = "_%d_%d" % (self.u, self.v)
        
        pos_suffix = lambda pos : "_%d_%d_%d" % (pos.x, pos.y, pos.i)
        cls_suffix = lambda cls : "_%d_%d" % (cls.x, cls.y)

        u_pos_var = lambda cp : "u_" + u + pos_suffix(cp.u)
        v_pos_var = lambda cp : "u_" + v + pos_suffix(cp.v)
        cp_var = lambda cp : 'e' + e + pos_suffix(cp.u) + pos_suffix(cp.v)

        u_cls_var = lambda stc : "u_" + u + cls_suffix(stc.u)
        v_cls_var = lambda stc : "u_" + v + cls_suffix(stc.v)
        stc_var = lambda stc : 'e' + e + cls_suffix(stc.u) + cls_suffix(stc.v)

        #Constraint declaration:

        self.cov_linearizations = []
        self.stretch_linearizations = []
        self.uncov_linearizations = []
        if DELAY_PIVOT:
            self.group_delays(fixed_nodes)
        else:
            self.td_assign = "- td" + e + " + "
            self.cov_assign = ""
        self.force_select = ""
        self.at_most_one_cov = ""
        self.uncov_assign = ""

        #Variable declaration:
        if not DELAY_PIVOT:
            self.bin_vars = set()
            #Otherwise, the set has already been declared.

        max_achievable_delay = -1
        min_achievable_delay = float('inf')
        for cp in self.cov_pairs:
            extend_var = cp_var(cp)
            #Linearizations:
            if self.u in fixed_nodes:
                extend_var = v_pos_var(cp)
            elif self.v in fixed_nodes:
                extend_var = u_pos_var(cp)
            else:
                self.cov_linearizations += linearize(u_pos_var(cp), v_pos_var(cp),\
                                                     cp_var(cp))
                self.bin_vars.add(u_pos_var(cp))
                self.bin_vars.add(v_pos_var(cp))
            self.bin_vars.add(extend_var) 
            
            if cp.td < min_achievable_delay:
                min_achievable_delay = cp.td
            if cp.td > max_achievable_delay:
                max_achievable_delay = cp.td

            if not DELAY_PIVOT:
                #Extend covered indicator assignment:
                #(We will simply take the uncovered indicator as the negation.)
                self.cov_assign += extend_var + " + "

                #Extend the delay assignment constraint:
                self.td_assign += print_float(cp.td) + ' ' + extend_var + " + "
        
        #Complete the covered indicator assignment:
        #(Also must bound the cover pair sum by 1)
        if not DELAY_PIVOT and self.cov_pairs:
            self.at_most_one_cov = self.cov_assign[:-len(" + ")]
            self.cov_assign = self.at_most_one_cov + " - cov_ind" + e + " = 0"
            self.at_most_one_cov += " <= 1"
            self.bin_vars.add("cov_ind" + e)

        #Now take a negation for the uncovered:
        if self.stretch_pairs and self.cov_pairs:
            self.uncov_assign = negate("cov_ind" + e, "uncov_ind" + e)
            self.bin_vars.add("cov_ind" + e)
            self.bin_vars.add("uncov_ind" + e)
            
        #Start the selection forcing constraint:
        if self.cov_pairs:
            self.force_select = "cov_ind" + e + " + "
        else:
            self.force_select = ""

        lin_prefix = 'u' if self.cov_pairs else ''
        #Split the uncov delay from the rest of the sum.
        #This allows for e.g., lazy constraining of stretches.
        if not DELAY_PIVOT and self.stretch_pairs:
            #FIXME: Maybe do the same for delay pivoting.
            if self.cov_pairs:
                self.td_assign += "elong" + e + " + " + print_float(self.init_td)\
                               + " uncov_ind" + e + " = 0"
            else:
                self.td_assign += "elong" + e + " = -" + print_float(self.init_td)
            self.elong_assign = "- elong" + e

        all_zero = True #Indicates that no stretch can bring a delay change.
        for stc in self.stretch_pairs:
            #Linearizations:
            self.stretch_linearizations += linearize(u_cls_var(stc),\
                                                     v_cls_var(stc),\
                                                     stc_var(stc))
            self.bin_vars.add(u_cls_var(stc))
            self.bin_vars.add(v_cls_var(stc))
            self.bin_vars.add(stc_var(stc))

            if self.cov_pairs:
                #Uncovered indicator product linearizations:
                self.uncov_linearizations += linearize(stc_var(stc),\
                                                       "uncov_ind" + e,\
                                                       'u' + stc_var(stc))
                self.bin_vars.add('u' + stc_var(stc))

            #Extend the selection forcing constraint:
            self.force_select += lin_prefix + stc_var(stc) + " + "

            if stc.td < min_achievable_delay:
                min_achievable_delay = stc.td
            if stc.td > max_achievable_delay:
                max_achievable_delay = stc.td

            if not DELAY_PIVOT:
                #Extend the delay assignment constraint:
                elong_inc = stc.td - self.init_td
                if all(c in ['0', '.'] for c in print_float(elong_inc)):
                    continue
                else:
                    all_zero = False
                if elong_inc < 0:
                    self.elong_assign += " - " + print_float(abs(elong_inc))\
                                       + ' ' + lin_prefix + stc_var(stc)
                else:
                    self.elong_assign += " + " + print_float(abs(elong_inc))\
                                       + ' ' + lin_prefix + stc_var(stc)

        #Complete the selection forcing constraint:
        self.force_select = self.force_select[:-len(" + ")] + " = 1"
    
        if not DELAY_PIVOT:
            #Complete the delay assignment constraint:
            if not self.stretch_pairs:
                self.td_assign = self.td_assign[:-len(" + ")] + " = 0"
            else:
                if all_zero:
                    del self.elong_assign
                else:
                    self.elong_assign += " = 0"

        #Generate the delay bound:
        if TIGHT_TD_BOUNDS:
            self.td_bound = print_float(min_achievable_delay) + " <= td" + e\
                          + " <= " + print_float(min(self.td_max, max_achievable_delay))
        else:
            self.td_bound = "td" + e + " <= " + print_float(self.td_max)
        if not DELAY_PIVOT and self.stretch_pairs:
            if all_zero:
                self.elong_bound = "elong" + e + " = 0"
            else:
                max_elong = max([stc.td for stc in self.stretch_pairs]) - self.init_td
                self.elong_bound = "elong" + e + " <= " + print_float(max_elong)

        self.td_bound_num = (min_achievable_delay, max_achievable_delay)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def print_lazy_constraints(self):
        """Prints all constraints that should be evaluated lazily.
 
        Parameters
        ----------
        None

        Returns
        -------
        str
            Text of the constraint block.
        """

        txt = "\* Lazy  (%s -> %s)   *\\\n" % (self.u_alias, self.v_alias)

        if hasattr(self, "elong_assign"):
            txt += self.elong_assign + "\n"
        #txt += self.force_select + "\n"
        #for lin in self.stretch_linearizations:
        #    txt += lin + "\n"
        #for lin in self.uncov_linearizations:
        #    txt += lin + "\n"

        return txt
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def print_constraints(self):
        """Prints all constraints for the given edge in
        a single block.

        Parameters
        ----------
        None

        Returns
        -------
        str
            Text of the constraint block.
        """

        txt = "\*   (%s -> %s)   *\\\n" % (self.u_alias, self.v_alias)
        txt += self.td_assign + "\n"
        if self.cov_pairs:
            #txt += self.at_most_one_cov + "\n"
            txt += self.cov_assign + "\n"
        if self.stretch_pairs and self.cov_pairs:
            txt += self.uncov_assign + "\n"
        for lin in self.cov_linearizations:
            txt += lin + "\n"
        
        txt += self.force_select + "\n"
        for lin in self.stretch_linearizations:
            txt += lin + "\n"
        for lin in self.uncov_linearizations:
            txt += lin + "\n"

        return txt
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def print_td_bound(self):
        """Prints the delay bound.

        Parameters
        ----------
        None
        
        Returns
        -------
        str
            Text of the bound.
        """

        txt = "\*   (%s -> %s)   *\\\n" % (self.u_alias, self.v_alias)
        txt += self.td_bound + "\n"
        if hasattr(self, "elong_bound"):
            txt += self.elong_bound + "\n" 

        return txt
    #-----------------------------------------------------------------------#
##########################################################################

##########################################################################
class Node(object):
    """Class describing a node of the timing graph which should
    enter the ILP formulation.

    Parameters
    ----------
    u : int
        Node identifier.

    Attributes
    ----------
    *Parameters
    positions : List[NodePos]
        List of all exact positions that the node can take.
    clusters : List[NodeCls]
        List of all cluster-level positions that the node can take.
    
    Methods
    -------
    add_position(x : int, y : int, i : int)
        Adds a new exact position.
    add_cluster(x : int, y: int)
        Adds a new cluster.
    set_alias(alias : str)
        Sets the alias of the node.
    set_init_cluster(cls : NodeCls)
        Sets the initial cluster of the node.
    generate_constraints()
        Generates constraints for the node.
    print_constraints()
        Prints all constraints for the given node in
        a single block.
    """

    #-----------------------------------------------------------------------#
    def __init__(self, u):
        """Constructor of the Node class.
        """

        self.u = u
        self.alias = str(u)
        self.positions = []
        self.clusters = []
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def add_position(self, x, y, i):
        """ Adds a new exact position.

        Also adds a new cluster if needed, as choosing a position
        implies choosing a cluster as well.

        Parameters
        ----------
        x : int
            x-coordinate of the cluster.
        y : int
            y-coordinate of the cluster.
        i : int
            LUT-index inside the cluster.

        Returns
        -------
            None
        """    

        position = NodePos(x = x, y = y, i = i)
        if not position in self.positions:
            self.positions.append(position)

        cluster = NodeCls(x = x, y = y)
        if not cluster in self.clusters:
            self.clusters.append(cluster)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def add_cluster(self, x, y):
        """ Adds a new cluster.

        Parameters
        ----------
        x : int
            x-coordinate of the cluster.
        y : int
            y-coordinate of the cluster.

        Returns
        -------
            None
        """

        cluster = NodeCls(x = x, y = y)
        if not cluster in self.clusters:
            self.clusters.append(cluster)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def set_alias(self, alias):
        """Sets the alias of the node.

        Parameters
        ----------
        alias : str

        Returns
        -------
        None
        """

        self.alias = alias
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def set_init_cluster(self, cls):
        """Sets the initial cluster of the node.
        
        Parameters
        ----------
        cls : NodeCls
            The initial cluster.

        Returns
        -------
        None
        """
    
        self.init_cls = cls
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def generate_constraints(self):
        """Generates constraints for the node.
        
        The following types of constraints are generated:

        1) Node must reside in exactly one cluster.
        2) Within one cluster, node can assume at most one
           exact position.
        3) Cluster membership whenever an exact position is chosen.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        #String conversions:

        pos_var = lambda pos : "u_%d_%d_%d_%d" % (self.u, pos.x, pos.y, pos.i)
        cls_var = lambda cls : "u_%d_%d_%d" % (self.u, cls.x, cls.y)

        #Constraint declaration:
        
        self.one_cluster = ""
        self.at_most_one_pos = {}
        self.cluster_membership = {}

        #Variable declaration:
        self.bin_vars = set()

        for pos in self.positions:
            cls = cls_var(pos)
            #At most one position per cluster:
            if cls in self.at_most_one_pos:
                self.at_most_one_pos[cls] += pos_var(pos) + " + "
            else:
                self.at_most_one_pos.update({cls : pos_var(pos) + " + "})
            #Position implies cluster membership:
            if cls in self.cluster_membership:
                self.cluster_membership[cls] += pos_var(pos) + " + "
            else:
                self.cluster_membership.update({cls : "- " + cls + " + "\
                                               + pos_var(pos) + " + "})
            self.bin_vars.add(pos_var(pos))
            self.bin_vars.add(cls)

        #Add a quadratic displacement tracking expression which can be used
        #as a secondary objective to help minimize the symmetries. No new
        #variables are needed as all the expressions will simply be summed
        #up and passed as objective in the end.
        self.displacement = ""
        get_displacement = lambda cls : (cls.x - self.init_cls.x) ** 2\
                                      + (cls.y - self.init_cls.y) ** 2

        init_included = False
        for cls in self.clusters:
            var = cls_var(cls)
            #Node belongs to exactly one cluster:
            self.one_cluster += var + " + "
            self.displacement += "%d %s + " % (get_displacement(cls), var)
            self.bin_vars.add(var)
        
            #Finish ``at most one position per cluster'':
            if var in self.at_most_one_pos:
                self.at_most_one_pos[var] =\
                self.at_most_one_pos[var][:-len(" + ")] + " <= 1"
 
            #Finish ``position implies cluster membership'':
            if var in self.cluster_membership:
                if not ONLY_COV_CLS or cls == self.init_cls: 
                    self.cluster_membership[var] =\
                    self.cluster_membership[var][:-len(" + ")] + " <= 0"
                    init_included = True
                else:
                    self.cluster_membership[var] =\
                    self.cluster_membership[var][:-len(" + ")] + " = 0"
                           
        #Finish ``exactly one cluster'': 
        self.one_cluster = self.one_cluster[:-len(" + ")] + " = 1"

        if not init_included:
            print("Warning: Initial cluster not included for node %s (%d)"\
                  % (self.alias, self.u))
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def print_lazy_constraints(self):
        """Prints all constraints that should be evaluated lazily.
 
        Parameters
        ----------
        None

        Returns
        -------
        str
            Text of the constraint block.
        """

        txt = "\* Lazy  %s   *\\\n" % self.alias
        #txt += self.one_cluster + "\n"
        #for cls in self.cluster_membership:
        #    #txt += self.at_most_one_pos[cls] + "\n"
        #    txt += self.cluster_membership[cls] + "\n"

        return txt
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def print_constraints(self):
        """Prints all constraints for the given node in
        a single block.

        Parameters
        ----------
        None

        Returns
        -------
        str
            Text of the constraint block.
        """

        txt = "\*   %s   *\\\n" % self.alias
        txt += self.one_cluster + "\n"
        for cls in self.cluster_membership:
            txt += self.at_most_one_pos[cls] + "\n"
            txt += self.cluster_membership[cls] + "\n"

        return txt
    #-----------------------------------------------------------------------#
##########################################################################

##########################################################################
class Cluster(object):
    """Class describing a cluster of the FPGA grid which should
    enter the ILP formulation.

    Parameters
    ----------
    x : int
        x-coordinate.
    y : int
        y-coordinate.

    Attributes
    ----------
    *Parameters    
    nodes : List[int]
        A list of nodes that can potentially be part of this cluster.
    fixed_inputs : List[int]
        A list of (initially) external cluster inputs required by
        the fixed nodes of the cluster.
    fixed_feedbacks : List[int]
        A list of (initially) internal inputs required by the fixed
        nodes of the cluster.
    
    Methods
    -------
    add_node(u : int)
        Adds a potential node to the cluster.
    add_fixed_input(i : int)
        Adds an input required by a fixed node.
    add_fixed_feedback(f : int)
        Adds a feedback required by a fixed node.
    generate_constraints(fanins : Dict[int, Set[int]],
                         fixed_nodes : Dict[int, NodeCls])
        Generates constraints for the cluster.
    print_constraints()
        Prints all constraints for the given cluster in
        a single block.
    print_input_indicator_bounds()
        Prints all bounds on input indicators.
    """

    #-----------------------------------------------------------------------#
    def __init__(self, x, y):
        """Constructor of the Cluster class.
        """

        cls = NodeCls(x, y)
        self.x = x
        self.y = y
        self.nodes = []
        self.fixed_inputs = []
        self.fixed_feedbacks = []
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def add_node(self, u):
        """Adds a potential node to the cluster.

        Parameters
        ----------
        u : int
            Node.

        Returns
        -------
        None
        """
        
        self.nodes.append(u)
    #-----------------------------------------------------------------------#
    
    #-----------------------------------------------------------------------#
    def add_fixed_input(self, i):
        """Adds an input required by a fixed node.

        Parameters
        ----------
        i : int
            Input.
        
        Returns
        -------
        None
        """
        
        if not i in self.fixed_inputs:
            self.fixed_inputs.append(i)
    #-----------------------------------------------------------------------#
 
    #-----------------------------------------------------------------------#
    def add_fixed_feedback(self, f):
        """Adds a feedback required by a fixed node.

        Parameters
        ----------
        f : int
            Feedback.
        
        Returns
        -------
        None
        """
        
        if not f in self.fixed_feedbacks:
            self.fixed_feedbacks.append(f)
    #-----------------------------------------------------------------------#
   
    #-----------------------------------------------------------------------#
    def generate_constraints(self, fanins, fixed_nodes):
        """Generates constraints for the cluster.

        The following types of constraints are generated:

        1) Bound on the maximum number of nodes in the cluster.
        2) Bound on the maximum number of inputs to the cluster.
        3) For each potential input entering 2), a binary indicator
           of presence of its children in the cluster.
        4) For each potential input entering 2), negation of its
           presence in the current cluster (if this is not a constant).
        5) Linearizations of the products of 3) and 4).

        Parameters
        ----------
        fanins : Dict[int, Set[int]]
            A dictionary of fanins for each node.
        fixed_nodes : Dict[int, NodeCls]
            A dictionary of cluster coordinates for each fixed node.
        
        Returns
        -------
        None
        """

        #NOTE: Even if there are no nodes in the cluster, we must generate
        #      constraints, as it could happen that some feedbacks will go
        #      away during moves.

        #.......................................................................#
        def group_children(self, fanins, fixed_nodes):
            """Goes through the set of potential cluster inputs
            and groups their children together.

            Parameters
            ----------
            fanins : Dict[int, Set[int]]
                A dictionary of fanins for each node.
            fixed_nodes : Dict[int, NodeCls]
                A dictionary of cluster coordinates for each fixed node.
            
            Returns
            -------
            Dict[int, Set[int]]
                A dictionary indexed by potential cluster inputs,
                holding sets of nodes adjacent to them in the current
                cluster.
            """

            adj = {}
            for v in self.nodes:
                for u in fanins.get(v, []):
                    if fixed_nodes.get(u, (-1, -1)) == (self.x, self.y):
                        continue
                    if u in self.fixed_inputs:
                        continue
                    try:
                        adj[u].add(v)
                    except:
                        adj.update({u : set([v])})

            return adj
        #.......................................................................#

        #String conversions:
        
        suffix = lambda u : "_%d_%d_%d" % (u, self.x, self.y)
        node_cls_var = lambda u : 'u' + suffix(u)
        input_fanout_var = lambda i : 'f' + suffix(i)
        input_fanout_bin_var = lambda i : "fb" + suffix(i)
        input_external_var = lambda i : "nu" + suffix(i)
        input_var = lambda i : 'i' + suffix(i)

        #Constraint declaration:

        self.max_inputs = I
        self.max_nodes = N - len([u for u in fixed_nodes\
                                  if fixed_nodes[u] == NodeCls(self.x, self.y)])
        
        self.node_limit = ""
        self.input_fanouts = {}
        self.input_fanouts_bin = {}
        self.input_presence_neg = {}
        self.input_linearizations = {}
        self.fixed_feedback_neg = {}
        self.input_limit = ""

        #Bound declaration:
        self.input_fanout_bounds = {}

        #Variable declaration:
        self.bin_vars = set()
        self.general_vars = set()

        #Node limit:
        for u in self.nodes:
            self.node_limit += node_cls_var(u) + " + "
            self.bin_vars.add(node_cls_var(u))
        if self.node_limit:
            self.node_limit = self.node_limit[:-len(" + ")]\
                            + " <= %d" % self.max_nodes
        
        fanouts = group_children(self, fanins, fixed_nodes)
        for i in sorted(fanouts):
            #Input fanout indication:
            self.input_fanouts.update({i : "- " + input_fanout_var(i) + " + "})
            self.general_vars.add(input_fanout_var(i))
            for u in sorted(fanouts[i]):
                self.input_fanouts[i] += node_cls_var(u) + " + "
                self.bin_vars.add(node_cls_var(u))
            self.input_fanouts[i] = self.input_fanouts[i][:-len(" + ")]\
                                    + " = 0"

            #Input fanout binarization:
            #fanout_lim = N
            fanout_lim = len(fanouts[i])
            #FIXME: In case fanout is limited to 1, we need no binarization nor
            #       a new indicator variable. Not crucial...
            self.input_fanout_bounds.update({i : fanout_lim})
            self.input_fanouts_bin.update({i : binarize(input_fanout_var(i),\
                                                        input_fanout_bin_var(i),\
                                                        fanout_lim)})
            self.bin_vars.add(input_fanout_bin_var(i))

            #Input not in cluster:
            if i in self.nodes:
                self.input_presence_neg.update({i : negate(node_cls_var(i),\
                                                           input_external_var(i))})
                self.bin_vars.add(input_external_var(i))

                #Input not in cluster AND its children in cluster:
                self.input_linearizations.update({i : linearize(input_external_var(i),\
                                                                input_fanout_bin_var(i),\
                                                                input_var(i))})
                self.bin_vars.add(input_var(i))
    
                #Extend the input count constraint:
                self.input_limit += input_var(i) + " + "
            else:
                #If the input can not enter the cluster, we can use the
                #binarized fanout indicator immediately.
                #We update the dictionaries with empty constraints to
                #maintain consistency.
                self.input_presence_neg.update({i : ""})
                self.input_linearizations.update({i : []})
                self.input_limit += input_fanout_bin_var(i) + " + "
        
        #Add feedbacks required by fixed nodes that can become external:    
        for f in self.fixed_feedbacks:
            if not f in self.nodes:
                #It may happen that in order to meet the timing requirements
                #for some connection, a feedback must be converted to an 
                #external input. Note that this is never the case if "f" is
                #fixed.
                continue
            self.fixed_feedback_neg.update({f : negate(node_cls_var(f),\
                                                       input_external_var(f))})
            self.bin_vars.add(node_cls_var(f))
            self.bin_vars.add(input_external_var(f))
            self.input_limit += input_external_var(f) + " + "
        if self.input_limit:
            self.input_limit = self.input_limit[:-len(" + ")] + " - "

        #Deduct external inputs required by fixed nodes
        #that can become internal:
        for i in self.fixed_inputs:
            if not i in self.nodes:
                #This input will stay external.
                self.max_inputs -= 1
                continue
            self.input_limit += node_cls_var(i) + " - "
            self.bin_vars.add(node_cls_var(i))
      
        if self.input_limit: 
            self.input_limit = self.input_limit[:-len(" - ")]\
                             + " <= %d" % self.max_inputs
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def print_constraints(self):
        """Prints all constraints for the given cluster in
        a single block.

        Parameters
        ----------
        None

        Returns
        -------
        str
            Text of the constraints.
        """

        txt = "\*   (%d, %d)   *\\\n" % (self.x, self.y)
        txt += self.node_limit + "\n"
        txt += self.input_limit + "\n"
        for i in self.input_fanouts:
            txt += self.input_fanouts[i] + "\n"
            for cstr in self.input_fanouts_bin[i]:
                txt += cstr + "\n"
            txt += self.input_presence_neg[i] + "\n"
            for lin in self.input_linearizations[i]:
                txt += lin + "\n"
        for f in self.fixed_feedback_neg:
            txt += self.fixed_feedback_neg[f] + "\n"

        return txt
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def print_input_indicator_bounds(self):
        """Prints all bounds on input indicators.

        Parameters
        ----------
        None

        Returns
        -------
        str
            Text of the bounds.
        """
        
        txt = "\*   (%d, %d)   *\\\n" % (self.x, self.y)
        for i in sorted(self.input_fanout_bounds):
            txt += "f_%d_%d_%d <= %d\n" % (i, self.x, self.y,\
                                           self.input_fanout_bounds[i])

        return txt
    #-----------------------------------------------------------------------#
##########################################################################

##########################################################################
class TimingGraph(object):
    """Class describing the timing graph.

    Parameters
    ----------
    tg : networkx.DiGraph
        Timing graph with each edge annotated by its delay.
    ch : Set[Tuple[str]]
        Set of edges whose delay will change.
    tarrs : Dict[str, float]
        A dictionary of arrival times, for those nodes whose arrival
        time can not change after ILP solving, but are parents of
        other nodes whose arrival time can change. All FF outputs
        are assigned a zero arrival time.
    node_counts : Dict[str, int]
        A mapping between node names and their counts used in Node
        and Edge classes. Here we need the strings, to be able to
        recognize the FFs. Note that multiple timing graph nodes
        map to the same number, because the numbers designate
        objects that are moved, which are in fact BLEs.

    Attributes
    ----------
    *Parameters

    Methods
    -------
    set_max_cpd(cpd : float)
        Sets the maximum arrival time.
    generate_node_constraints(u : str)
        Generates constraints for a single node.
    generate_constraints()
        Generates all constraints.
    print_constraints()
        Prints all constraints.
    print_bounds()
        Prints all bounds.
    """

    #-----------------------------------------------------------------------#
    def __init__(self, tg, ch, tarrs, node_counts):
        """Constructor of the TimingGraph class.
        """

        self.tg = tg
        self.ch = ch
        self.tarrs = tarrs
        self.node_counts = node_counts

        self.node_constraints = {}
        self.node_bounds = {}
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def set_max_cpd(self, cpd):
        """Sets the maximum arrival time.

        Parameters
        ----------
        cpd : float

        Returns
        -------
        None
        """

        self.cpd = cpd
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def generate_node_constraints(self, u):
        """Generates node constraints for a single node.

        Parameters
        ----------
        u : str
            Name of the node.
        
        Returns
        -------
        None
        """

        u_tarr = self.tarrs.get(u, -1)
        if u_tarr >= 0:
            return

        #String conversions:
        node_var = lambda u : "ta_%d" % self.node_counts[u]
        edge_var = lambda u, v : "td_%d_%d" % (self.node_counts[u],\
                                               self.node_counts[v])
        csts = []
        lb = -1
        indeg = self.tg.in_degree(u)
        rel = " = " if indeg == 1 else " >= "
        for p in self.tg.pred[u]:
            if p in self.tarrs:
                if not (p, u) in self.ch:
                    td = self.tarrs[p] + self.tg[p][u]["td"]
                    lb = max(td, lb)
                    continue
                cst = node_var(u) + " - " + edge_var(p, u)\
                    + rel + print_float(self.tarrs[p])
            elif not (p, u) in self.ch:
                cst = node_var(u) + " - " + node_var(p)\
                    + rel + print_float(self.tg[p][u]["td"])
            else:
               cst = node_var(u) + " - " + node_var(p) + " - "\
                   + edge_var(p, u) + rel + '0'
            csts.append(cst)

        self.node_constraints.update({u : csts})
        if lb >= 0:
            self.node_bounds.update({u : lb})
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def remove_irrelevant_edges(self, max_edge_delays):
        """Removes the edges that can not become critical even if all
        connections are as slow as possible.

        Parameters
        ----------
        max_edge_delays : Dict[Tuple[str], float]
            Maximum delay of each changeable edge.

        Returns
        -------
        None
        """

        tg_ub = self.tg.copy()
        for e in max_edge_delays:
            u, v = e
            tg_ub[u][v]["td"] = max_edge_delays[e]
        
        topo_nodes = list(nx.topological_sort(tg_ub))
        rev_topo_nodes = reversed(topo_nodes)
    
        for v in topo_nodes:
            tg_ub.node[v]["ta"] = max([tg_ub.node[u]["ta"] + tg_ub[u][v]["td"]\
                                    for u in tg_ub.pred[v]] + [0])
        for u in rev_topo_nodes:
            tg_ub.node[u]["tr"] = min([tg_ub.node[v]["tr"] - tg_ub[u][v]["td"]\
                                    for v in tg_ub[u]] + [self.pruning_cpd])

        total = 0
        irrelevant = []
        for u, v in tg_ub.edges():
            slack = tg_ub.node[v]["tr"] - tg_ub.node[u]["ta"] - tg_ub[u][v]["td"]
            total += 1
            if slack >= 0:
                irrelevant.append((u, v))

        print "All edges: ", self.tg.number_of_edges()
        for u, v in irrelevant:
            self.tg.remove_edge(u, v)
        isolated_nodes = [u for u in self.tg if self.tg.in_degree(u) == 0\
                          and self.tg.out_degree(u) == 0]
        for u in isolated_nodes:
            self.tg.remove_node(u)
        print "Left edges: ", self.tg.number_of_edges()

        self.generate_constraints()
    #-----------------------------------------------------------------------#
    
    #-----------------------------------------------------------------------#
    def tighten_arrival_bounds(self, min_edge_delays, max_edge_delays):
        """Performs arrival time bound tightening, by doing STA with all
        connections being as fast as possible (minimum arrival times) and
        all connections being as slow as possible (maximum arrival times).

        Parameters
        ----------
        min_edge_delays : Dict[Tuple[str], float]
            Minimum delay of each changeable edge.
        max_edge_delays : Dict[Tuple[str], float]
            Maximum delay of each changeable edge.

        Returns
        -------
        bool
            True if the bounds are consistent with the target CPD.
            False if there is negative slack under least possible edge delays.
        """
        
        tg_lb = self.tg.copy()
        tg_ub = self.tg.copy()

        for e in min_edge_delays:
            u, v = e
            if not tg_lb.has_edge(u, v):
                print(u, v, "swept.")
                continue
            tg_lb[u][v]["td"] = min_edge_delays[e]
            tg_ub[u][v]["td"] = max_edge_delays[e]

        #Round the delays so that accumulated error does not cause
        #infeasibility. WARNING: Be careful about that!
        for u, v in tg_lb.edges():
            tg_lb[u][v]["td"] = float(print_float(tg_lb[u][v]["td"]))
            tg_ub[u][v]["td"] = float(print_float(tg_ub[u][v]["td"]))

        topo_nodes = list(nx.topological_sort(tg_lb))
        rev_topo_nodes = reversed(topo_nodes)

        bounds = {}
        for v in topo_nodes:
            ta_lb = max([tg_lb.node[u]["ta"] + tg_lb[u][v]["td"]\
                         for u in tg_lb.pred[v]] + [0])
            ta_ub = max([tg_ub.node[u]["ta"] + tg_ub[u][v]["td"]\
                         for u in tg_ub.pred[v]] + [0])
            tg_lb.node[v]["ta"] = ta_lb
            tg_ub.node[v]["ta"] = ta_ub
            if not v in self.tarrs:
                bounds.update({v : (ta_lb, ta_ub)})
        for u in rev_topo_nodes:
            tr = min([tg_lb.node[v]["tr"] - tg_lb[u][v]["td"]\
                      for v in tg_lb[u]] + [self.cpd])
            if tr < 0:
                #Negative slack.
                return False

            tg_lb.node[u]["tr"] = tr
            ta_ub = bounds.get(u, (-1, -1))[1]
            if tr < ta_ub:
                bounds[u] = (bounds[u][0], tr)

        self.tight_node_bounds = bounds
        self.lb_cpd = max([tg_lb.node[u]["ta"] for u in topo_nodes])

        #Store the graphs for future use.
        self.tg_lb = tg_lb
        self.tg_ub = tg_ub

        return True
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

        self.node_constraints = {}
        self.node_bounds = {}
        for u in self.tg:
            self.generate_node_constraints(u)
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

        #String conversions:
        node_var = lambda u : "ta_%d" % self.node_counts[u]

        txt = ""
        for u in self.node_constraints:
            txt += "\*   %s   *\\\n" % u
            for cst in self.node_constraints[u]:
                txt += cst + "\n"
            if not self.tg.out_degree(u):
                txt += "Tamax - " + node_var(u) + " >= 0\n"

        return txt
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def print_bounds(self, forced_cpd_lb = None):
        """Prints all bounds.

        Parameters
        ----------
        forced_cpd_lb : Optional[float]
            An optional forced CPD lower bound.
            This is useful when solving a sequence of related ILPs.

        Returns
        -------
        str
            Text of the bounds.
        """

        #String conversions:
        node_var = lambda u : "ta_%d" % self.node_counts[u]

        txt = ""
        if hasattr(self, "tight_node_bounds"):
            for u in self.tight_node_bounds:
                txt += print_float(self.tight_node_bounds[u][0]) + " <= " + node_var(u)\
                     + " <= " + print_float(self.tight_node_bounds[u][1]) + "\n"
            lb_str = print_float(self.lb_cpd)
            if forced_cpd_lb is not None:
                lb_str = print_float(forced_cpd_lb)
            txt += lb_str + " <= Tamax <= " + print_float(self.cpd) + "\n"
        else:
            for u in self.node_bounds:
                txt += node_var(u) + " >= " + print_float(self.node_bounds[u]) + "\n"
            lb_str = print_float(max(list(self.tarrs.values()) + [0]))
            if forced_cpd_lb is not None:
                lb_str = print_float(forced_cpd_lb)
            txt += lb_str + " <= Tamax <= " + print_float(self.cpd) + "\n"

        return txt
    #-----------------------------------------------------------------------#
##########################################################################

##########################################################################
class ILP(object):
    """The main class describing the ILP generator.

    Parameters
    ----------
    tg : networkx.DiGraph
        Timing graph with each edge annotated by its delay.
    tarrs : Dict[str, float]
        A dictionary of arrival times, for those nodes whose arrival
        time can not change after ILP solving, but are parents of
        other nodes whose arrival time can change. All FF outputs
        are assigned a zero arrival time.
    node_counts : Dict[str, int]
        A mapping between node names and their counts used in Node
        and Edge classes. Here we need the strings, to be able to
        recognize the FFs. Note that multiple timing graph nodes
        map to the same number, because the numbers designate
        objects that are moved, which are in fact BLEs.
    init_positions : Dict[int, NodeCls]
        Initial clusters of the movable nodes.
        Necessary for cluster input counting.
    fixed_nodes : Dict[int, NodeCls]
        A mapping of nodes not movable by the ILP to their
        respective cluster positions. 
    cov_map : Dict[BasicEdge, List[CovPair]]
        A partial mapping of the timing-graph edge set to direct
        connections of the FPGA grid. Contains those edges whose
        delay is affected by the ILP and that can be covered by
        a direct connection.
    stretch_map : Dict[BasicEdge, List[StretchPair]]
        A partial mapping of the timing-graph edge set to cluster
        pairs of the FPGA grid. Contains those edges whose delay is
        affected by the ILP.
    max_tds : Dict[BasicEdge, float]
        Maximum tollerable delays for those edges that are affected
        by the ILP.
    maad : Dict[BasicEdge, float]
        Maximum achievable delays for those edges that are affected
        by the ILP.
    cpd : float
        The desired critical path delay.
    forced_cpd_lb : Optional[float]
        Optional CPD lower bound. This is useful when solving
        a sequence of related ILPs.
    opt_mode_override : Optional[str]
        Overrides the OPT_MODE parameter.
    pruning_cpd : float
        Minimum expected CPD which is to be used for pruning irrelevant
        edges. If the target CPD is used for pruning, any objective
        value lower than the target may be invalid, because some edges
        that may invalidate it could have been removed.
    read_mst : Optional[bool], default = False
        Decides if an attempt to read a starting solution is to be made.

    Attributes
    ----------
    *Parameters
    edges : List[Edge]
        A list of edge objects, for all entries in the union of cov_map
        and stretch_map.
    nodes : List[Node]
        A list of node objects, for all nodes touched by the edges in
        cov_map and stretch_map.
    clusters : List[Cluster]
        A list of cluster objects, for all clusters in which nodes
        from nodes can reside.
    
    Methods
    -------
    build_rev_adj()
        Converts the timing graph to reverse adjacency as
        required by other methods.
    build_ch()
        Constructs the set of changeable connections from the
        cov_map and stretch_map keys.
    build_edges()
        Constructs all edge objects.
    build_nodes()
        Construct all node objects.
    build_clusters()
        Constructs all cluster objects.
    build_tg()
        Constructs the timing graph object.
    build_all()
        Constructs all objects.
    set_incumbent_coverage(List[Tuple[str]])
        Sets the coverage of the incumbent solution that is needed if 
        maximizing coverage similarity is the objective.
    find_node(u : int or str)
        Finds a node object based on an identifier.
    find_edge(u : int or str, v : int or str)
        Finds an edge object based on endpoint identifiers.
    find_fanin(u : int or str)
        Returns incoming edge objects incident to u.
    find_fanout(u : int or str)
        Returns outgoing edge objects incident to u.
    alias_nodes(aliases : Dict[int, str])
        Assigns aliases to a subset of nodes and relabels
        edges accordnigly.
    generate_pos_overlap_constraints()
        Generates the exact position overlap removal constraints.
    generate_constraints()
        Generates constraints for all objects.
    print_constraints()
        Prints all constraints.
    print_bounds()
        Prints all bounds.
    print_variables()
        Prints all variables.
    print_ilp()
        Prints the entire ILP formulation.
    solve(timeout : int)
        Calls the appropriate solver to solve the ILP.
    """

    #-----------------------------------------------------------------------#
    def __init__(self, tg, tarrs, node_counts, init_positions, fixed_nodes,\
                 cov_map, stretch_map, max_tds, maad, cpd, forced_cpd_lb = None,\
                 opt_mode_override = None, pruning_cpd = None, read_mst = False):
        """Constructor of the ILP class.
        """

        if opt_mode_override is not None:
            global OPT_MODE
            OPT_MODE = opt_mode_override

        if pruning_cpd is None:
            self.pruning_cpd = cpd
        else:
            self.pruning_cpd = pruning_cpd

        self.read_mst = read_mst

        self.edges = []
        self.nodes = []
        self.clusters = []

        self.bin_vars = set()
        self.general_vars = set()

        self.init_positions = init_positions
        self.fixed_nodes = fixed_nodes
        self.cov_map = cov_map
        self.stretch_map = stretch_map
        self.max_tds = max_tds
        self.maad = maad

        #These parameters are checked by the TimingGraph constructor.
        self.tg = tg.copy()
        self.tarrs = tarrs
        self.node_counts = node_counts

        self.inv_node_counts = {self.node_counts[u] : u[:-1] + 'd'\
                               if u.endswith(".q") else u\
                               for u in self.node_counts}
        self.cpd = cpd
        self.forced_cpd_lb = forced_cpd_lb
        self.build_rev_adj()

        self.build_ch()
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def build_rev_adj(self):
        """Converts the timing graph to reverse adjacency as
        required by other methods.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.rev_adj = {self.node_counts[u] : [] for u in self.tg}
        for u, v in self.tg.edges():
            u_int = self.node_counts[u]
            v_int = self.node_counts[v]
            if u_int != v_int:
                #Loops are not influenced by moves, but including them
                #in the formulation would create problems for input counting.
                self.rev_adj[v_int].append(u_int)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def build_ch(self):
        """Constructs the set of changeable connections from the
        cov_map and stretch_map keys.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.ch = set()
        for e in self.cov_map:
            u_str = self.inv_node_counts[e.u]
            if u_str.endswith(".d"):
                u_str = u_str[:-1] + 'q'
            v_str = self.inv_node_counts[e.v]
            if v_str.endswith(".q"):
                v_str = v_str[:-1] + 'd'
            self.ch.add((u_str, v_str))
        for e in self.stretch_map:
            u_str = self.inv_node_counts[e.u]
            if u_str.endswith(".d"):
                u_str = u_str[:-1] + 'q'
            v_str = self.inv_node_counts[e.v]
            if v_str.endswith(".q"):
                v_str = v_str[:-1] + 'd'
            self.ch.add((u_str, v_str))
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def build_edges(self):
        """Constructs all edge objects.
        
        Parameters
        ----------
        None

        Returns
        -------
        Dict[int, List[NodePos]]
            A dictionary of exact node positions for a subset of nodes.
        Dict[int, List[NodeCls]]
            A dictionary of cluster node positions for a subset of nodes.
            Only those clusters for which exact positions were not
            already specified are included.
        """

        node_pos_dict = {}
        node_cls_dict = {}
        created_edges = []
        for e in self.cov_map:
            u, v = e
            td_max = self.max_tds[e]
            u_str = self.inv_node_counts[u]
            if u_str.endswith(".d"):
                u_str = u_str[:-1] + 'q'
            v_str = self.inv_node_counts[v]
            if v_str.endswith(".q"):
                v_str = v_str[:-1] + 'd'
            e_obj = Edge(u, v, self.tg[u_str][v_str]["td"])
            e_obj.set_max_td(td_max)
            for cov in self.cov_map[e]:
                if FILTER_SLOW and cov.td > td_max:
                    continue
                e_obj.add_cov(cov.u, cov.v, cov.td)
                if not u in self.fixed_nodes:
                    try:
                        if not cov.u in node_pos_dict[u]:
                            node_pos_dict[u].append(cov.u)
                    except:
                        node_pos_dict.update({u : [cov.u]})
                if not v in self.fixed_nodes:
                    try:
                        if not cov.v in node_pos_dict[v]:
                            node_pos_dict[v].append(cov.v)
                    except:
                        node_pos_dict.update({v : [cov.v]})
            if e_obj.cov_pairs:
                self.edges.append(e_obj)
                created_edges.append(e)

        for e in self.stretch_map:
            u, v = e
            new_edge = False
            td_max = self.max_tds[e]
            try:
                e_obj = self.edges[created_edges.index(e)]
            except:
                u_str = self.inv_node_counts[u]
                if u_str.endswith(".d"):
                    u_str = u_str[:-1] + 'q'
                v_str = self.inv_node_counts[v]
                if v_str.endswith(".q"):
                    v_str = v_str[:-1] + 'd'
                e_obj = Edge(u, v, self.tg[u_str][v_str]["td"])
                e_obj.set_max_td(td_max)
                new_edge = True
            for stretch in self.stretch_map[e]:
                if FILTER_SLOW and stretch.td > td_max:
                    continue
                e_obj.add_stretch(stretch.u, stretch.v, stretch.td)
                if not u in self.fixed_nodes:
                    try:
                        found = False
                        for pos in node_pos_dict.get(u, []):
                            if stretch.u == node_pos_to_node_cls(pos):
                                found = True
                                break
                        if not found and not stretch.u in node_cls_dict[u]:
                            node_cls_dict[u].append(stretch.u)
                    except:
                        found = False
                        for pos in node_pos_dict.get(u, []):
                            if stretch.u == node_pos_to_node_cls(pos):
                                found = True
                                break
                        if not found:
                            node_cls_dict.update({u : [stretch.u]})
                if not v in self.fixed_nodes:
                    try:
                        found = False
                        for pos in node_pos_dict.get(v, []):
                            if stretch.v == node_pos_to_node_cls(pos):
                                found = True
                                break
                        if not found and not stretch.v in node_cls_dict[v]:
                            node_cls_dict[v].append(stretch.v)
                    except:
                        found = False
                        for pos in node_pos_dict.get(v, []):
                            if stretch.v == node_pos_to_node_cls(pos):
                                found = True
                                break
                        if not found:
                            node_cls_dict.update({v : [stretch.v]})
            if e_obj.stretch_pairs and new_edge:    
                self.edges.append(e_obj)
                created_edges.append(e)

        #Keeping the lists sorted helps potential debug.
        for u in node_pos_dict:
            node_pos_dict[u].sort()
        for u in node_cls_dict:
            node_cls_dict[u].sort()

        return node_pos_dict, node_cls_dict
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def build_nodes(self, node_pos_dict, node_cls_dict, init_clusters):
        """Constructs all node objects.
        
        Parameters
        ----------
        node_pos_dict : Dict[int, List[NodePos]]
            A dictionary of exact node positions for a subset of nodes.
        node_cls_dict : Dict[int, List[NodeCls]]
            A dictionary of cluster node positions for a subset of nodes.
            Only those clusters for which exact positions were not
            already specified are included.
        init_clusters : Dict[int, NodeCls]
            A dictionary of initial clusters for all nodes.

        Returns
        -------
        Dict[NodeCls, List[int]]
            A dictionary of nodes belonging to each cluster, indexed
            by the cluster's coordinates.
        """

        cls_dict = {}
        pos_keys = sorted(node_pos_dict)
        for u in pos_keys:
            u_obj = Node(u)
            u_obj.set_init_cluster(init_clusters[u])
            for pos in node_pos_dict[u]:
                u_obj.add_position(pos.x, pos.y, pos.i)
                cls = NodeCls(pos.x, pos.y)
                try:
                    if not u in cls_dict[cls]:
                        cls_dict[cls].append(u)
                except:
                    cls_dict.update({cls : [u]})
            self.nodes.append(u_obj)

        cls_keys = sorted(node_cls_dict)
        for u in cls_keys:
            new_node = False
            try:
                u_obj = self.nodes[pos_keys.index(u)]
            except:
                u_obj = Node(u)
                u_obj.set_init_cluster(init_clusters[u])
                new_node = True
            for cls in node_cls_dict[u]:
                u_obj.add_cluster(cls.x, cls.y)
                try:
                    #NOTE: Clusters are passed for only those
                    #      (node, cluster) pairs for which there
                    #      is no exact position. Hence, it can not
                    #      happen that we insert the same node twice.
                    cls_dict[cls].append(u)
                except:
                    cls_dict.update({cls : [u]})
            if new_node:
                self.nodes.append(u_obj)

        return cls_dict
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def count_fanouts(self):
        """Generates constraints that assign fanout counters for each node,
        counting the number of outgoing edges selected to be covered by the
        solver. This is useful for imposing position constraints.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        fanout_dict = {}
        for e in self.cov_map:
            try:
                fanout_dict[e.u].append(e.v)
            except:
                fanout_dict.update({e.u : [e.v]})

        #String conversions:

        cov_ind = lambda u, v : "cov_ind_%d_%d" % (u, v)
        fo_var = lambda u : "fo_%d" % u

        fo_assigns = {}
        fo_bounds = {} 
        for u in fanout_dict:
            fo_assign = "- " + fo_var(u) + " + "
            self.general_vars.add(fo_var(u))
            for v in fanout_dict[u]:
                fo_assign += cov_ind(u, v) + " + "
            fo_assigns.update({u : fo_assign[:-len(" + ")] + " = 0"})
            fo_bounds.update({u : "0 <= " + fo_var(u) + " <= %d" % len(fanout_dict[u])})

        self.fo_assigns = fo_assigns
        self.fo_bounds = fo_bounds
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def count_fanins(self):
        """Generates constraints that assign fanin counters for each node,
        counting the number of incoming edges selected to be covered by the
        solver. This is useful for imposing position constraints.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        fanin_dict = {}
        for e in self.cov_map:
            try:
                fanin_dict[e.v].append(e.u)
            except:
                fanin_dict.update({e.v : [e.u]})

        #String conversions:

        cov_ind = lambda u, v : "cov_ind_%d_%d" % (u, v)
        fi_var = lambda u : "fi_%d" % u

        fi_assigns = {}
        fi_bounds = {}
        for v in fanin_dict:
            fi_assign = "- " + fi_var(v) + " + "
            self.general_vars.add(fi_var(v))
            for u in fanin_dict[v]:
                fi_assign += cov_ind(u, v) + " + "
            fi_assigns.update({v : fi_assign[:-len(" + ")] + " = 0"})
            fi_bounds.update({v : "0 <= " + fi_var(v) + " <= %d" % len(fanin_dict[v])})

        self.fi_assigns = fi_assigns
        self.fi_bounds = fi_bounds
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def infer_pattern_nodes(self):
        """Generates constraints that extract pattern node for each
        movable node, based on its position variable.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        #String conversions:
        
        pos_var = lambda u, pos : "u_%d_%d_%d_%d" % (u, pos.x, pos.y, pos.i)
        pat_pos_var = lambda u, pos : "u_%d_%d" % (u, pos.i)

        pat_pos_assigns = {}
        for u_obj in self.nodes:
            u = u_obj.u
            pos_dict = {}
            pat_pos_assigns.update({u : ""})
            for pos in u_obj.positions:
                pat_pos = pat_pos_var(u, pos)
                self.bin_vars.add(pat_pos)
                try:
                    pos_dict[pat_pos].append(pos)
                except:
                    pos_dict.update({pat_pos : [pos]})
            for pat_pos in pos_dict:
                pat_pos_assign = "- " + pat_pos + " + "
                for pos in pos_dict[pat_pos]:
                    pat_pos_assign += pos_var(u, pos) + " + "
                pat_pos_assigns[u] += pat_pos_assign[:-len(" + ")] + " = 0\n"
    
        self.pat_pos_assigns = pat_pos_assigns 
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def load_pattern_edges(self, pattern_edges):
        """Loads the pattern edges.

        Parameters
        ----------
        pattern_edges : List[Tuple[int]]
            A list of all pattern edges, without offsets (only LUT indices).
        
        Returns
        -------
        None
        """
        
        self.pattern_edges = pattern_edges
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def fanout_implications(self):
        """Generates constraints on pattern position variables, depending
        on the current value of the fanout counter.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        fanout_dict = {}
        for u, v in self.pattern_edges:
            try:
                fanout_dict[u] += 1
            except:
                fanout_dict[u] = 1
            if not v in fanout_dict:
                fanout_dict.update({v : 0})

        fanout_dict = {fo : [u for u in fanout_dict if fanout_dict[u] == fo]\
                       for fo in fanout_dict.values()}
        available_fanouts = sorted(fanout_dict.keys(), reverse = True)

        #String conversions:

        fo_var = lambda u : "fo_%d" % u
        fo_bvar = lambda u, val : "fob_%d_%d" % (u, val)
        pos_var = lambda u, p : "u_%d_%d" % (u, p)

        fanout_implications = {}
        rem_list = []
        for u in self.fo_bounds:
            lines = self.pat_pos_assigns[u].splitlines()
            positions = [int(l.split(" + ")[0].split("- u_%d_" % u)[1])\
                         for l in lines]
            csts = ""
            max_fo = int(self.fo_bounds[u].rsplit(" <= ", 1)[1])
            if max_fo < 2:
                #NOTE: If fanout is smaller than 2, the constraints are
                #completely reduntant, as there would be no cov pair to
                #force the node to the banned positions in any case.
                rem_list.append(u)
                continue
            for fo in range(2, max_fo + 1):
                valid_positions = []
                banned_exist = False
                for cmp_fo in available_fanouts:
                    if cmp_fo < fo:
                        banned_exist = True
                        break
                    for pos in fanout_dict[cmp_fo]:
                        if pos in positions:
                            valid_positions.append(pos)
                if not banned_exist:
                    continue
                if not valid_positions:
                    prev_bound = self.fo_bounds[u]
                    new_bound = prev_bound.rsplit(" <= ", 1)[0] + " <= %d"\
                              % (fo - 1)
                    self.fo_bounds[u] = new_bound
                    break
                if max_fo > 1:
                    geq_csts = geq_indicate(fo_var(u), fo_bvar(u, fo), max_fo, fo)
                    self.bin_vars.add(fo_bvar(u, fo))
                    cst = geq_csts[0] + "\n" + geq_csts[1] + "\n"
                    cst += "- " + fo_bvar(u, fo) + " + "
                else:
                    cst = "- " + fo_var(u) + " + "
                for pos in valid_positions:
                    cst += pos_var(u, pos) + " + "
                cst = cst[:-len(" + ")] + " >= 0\n"
                csts += cst
            if not csts:
                rem_list.append(u)
                continue
            fanout_implications.update({u : csts})

        for u in rem_list:
            del self.fo_assigns[u]
            del self.fo_bounds[u]
            self.general_vars.remove("fo_%d" % u)

        self.fanout_implications = fanout_implications
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def fanin_implications(self):
        """Generates constraints on pattern position variables, depending
        on the current value of the fanin counter.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        fanin_dict = {}
        for u, v in self.pattern_edges:
            try:
                fanin_dict[v] += 1
            except:
                fanin_dict[v] = 1
            if not u in fanin_dict:
                fanin_dict.update({u : 0})

        fanin_dict = {fi : [u for u in fanin_dict if fanin_dict[u] == fi]\
                       for fi in fanin_dict.values()}
        available_fanins = sorted(fanin_dict.keys(), reverse = True)

        #String conversions:

        fi_var = lambda u : "fi_%d" % u
        fi_bvar = lambda u, val : "fib_%d_%d" % (u, val)
        pos_var = lambda u, p : "u_%d_%d" % (u, p)

        fanin_implications = {}
        rem_list = []
        for u in self.fi_bounds:
            lines = self.pat_pos_assigns[u].splitlines()
            positions = [int(l.split(" + ")[0].split("- u_%d_" % u)[1])\
                         for l in lines]
            csts = ""
            max_fi = int(self.fi_bounds[u].rsplit(" <= ", 1)[1])
            if max_fi < 2:
                #NOTE: If fanin is smaller than 2, the constraints are
                #completely reduntant, as there would be no cov pair to
                #force the node to the banned positions in any case.
                rem_list.append(u)
                continue
            for fi in range(2, max_fi + 1):
                valid_positions = []
                banned_exist = False
                for cmp_fi in available_fanins:
                    if cmp_fi < fi:
                        banned_exist = True
                        break
                    for pos in fanin_dict[cmp_fi]:
                        if pos in positions:
                            valid_positions.append(pos)
                if not banned_exist:
                    continue
                if not valid_positions:
                    prev_bound = self.fi_bounds[u]
                    new_bound = prev_bound.rsplit(" <= ", 1)[0] + " <= %d"\
                              % (fi - 1)
                    self.fi_bounds[u] = new_bound
                    break
                if max_fi > 1:
                    geq_csts = geq_indicate(fi_var(u), fi_bvar(u, fi), max_fi, fi)
                    self.bin_vars.add(fi_bvar(u, fi))
                    cst = geq_csts[0] + "\n" + geq_csts[1] + "\n"
                    cst += "- " + fi_bvar(u, fi) + " + "
                else:
                    cst = "- " + fi_var(u) + " + "
                for pos in valid_positions:
                    cst += pos_var(u, pos) + " + "
                cst = cst[:-len(" + ")] + " >= 0\n"
                csts += cst
            if not csts:
                rem_list.append(u)
                continue
            fanin_implications.update({u : csts})

        for u in rem_list:
            del self.fi_assigns[u]
            del self.fi_bounds[u]
            self.general_vars.remove("fi_%d" % u)

        self.fanin_implications = fanin_implications
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def build_clusters(self, cls_dict):
        """Constructs all cluster objects.
        
        Parameters
        ----------
        Dict[NodeCls, List[int]]
            A dictionary of nodes belonging to each cluster, indexed
            by the cluster's coordinates.

        Returns
        -------
        None
        """

        for cls in cls_dict:
            cls_obj = Cluster(cls.x, cls.y)
            for u in cls_dict[cls]:
                cls_obj.add_node(u)
            for v in self.fixed_nodes:
                if self.fixed_nodes[v] != cls:
                    continue
                for u in self.rev_adj[v]:
                    if self.fixed_nodes.get(u, None) == cls:
                        continue
                    u_cls = self.fixed_nodes.get(u, None)
                    if not u_cls is None:
                        cls_obj.add_fixed_input(u)
                        continue
                    if self.init_positions[u] == cls:
                        cls_obj.add_fixed_feedback(u)
                    else:
                        cls_obj.add_fixed_input(u)
            self.clusters.append(cls_obj)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def build_tg(self):
        """Constructs the timing graph object.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.timing = TimingGraph(self.tg, self.ch, self.tarrs, self.node_counts)
        self.timing.set_max_cpd(self.cpd)

        self.timing.pruning_cpd = self.pruning_cpd
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def build_all(self):
        """Constructs all objects.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
    
        node_pos_dict, node_cls_dict = self.build_edges()
        init_clusters = {self.node_counts[u] : self.tg.node[u]["coords"]\
                         for u in self.tg}
        cls_dict = self.build_nodes(node_pos_dict, node_cls_dict, init_clusters)
        if ENFORCE_LEGALIZATION:
            self.build_clusters(cls_dict)
        self.build_tg()
    #-----------------------------------------------------------------------#
    
    #-----------------------------------------------------------------------#
    def set_incumbent_coverage(self, coverage):
        """Sets the coverage of the incumbent solution that is needed
        if maximizing coverage similarity is the objective.

        Parameters
        ----------
        coverage : List[Tuple[str]]
            A list of covered edges.

        Returns
        -------
        None
        
        Raises
        ------
        ValueError
            If the Edge objects have not been built yet.
        """

        if not self.edges:
            print("Edge objects must be built before setting the incumbent"\
                  + " coverage, to ensure that all coverage indices are"\
                  + " declared as binary.")
            raise ValueError

        #String conversion:
        e_var = lambda e : "cov_ind_%d_%d" % (self.node_counts.get(e[0], -1),\
                                              self.node_counts.get(e[1], -1))

        available = [var for e in self.edges for var in e.bin_vars\
                     if var.startswith("cov_ind")]
        
        cst = ""
        for e in coverage:
            ind = e_var(e)
            if ind in available:
                cst += ind + " + "
        
        self.incumbent_coverage = cst[:-len(" + ")] if cst else "NOTHING"
        
        #Set minimum coverage.
        if FORCE_COVERAGE_INCREASE:
            cst = ""
            for cov_ind in available:
                cst += cov_ind + " + "
            cst = cst[:-len(" + ")] + " >= %d\n" % len(coverage)
            self.least_coverage = cst        
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def find_node(self, u):
        """Finds a node object based on an identifier.

        Parameters
        ----------
        u : int or str
            Node identifier.
            If a string is passed, an alias-based lookup is performed.
    
        Returns
        -------
        Node or None
            The found object, or None if nothing is found.
        """

        for u_obj in self.nodes:
            if isinstance(u, str):
                if u == u_obj.alias:
                    return u_obj
            else:
                if u == u_obj.u:
                    return u_obj
            
        return None
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def find_edge(self, u, v):
        """Finds an edge object based on endpoint identifiers.

        Parameters
        ----------
        u : int or str
            Tail node identifier.
            If a string is passed, an alias-based lookup is performed.
 
        v : int or str
            Head node identifier.
            If a string is passed, an alias-based lookup is performed.
   
        Returns
        -------
        Node or None
            The found object, or None if nothing is found.
        """

        for e_obj in self.edges:
            if isinstance(u, str):
                if u == e_obj.u_alias:
                    if isinstance(v, str):
                        if v == e_obj.v_alias:
                            return e_obj
                    else:
                        if v == e_obj.v:
                            return e_obj
            else:
                if u == e_obj.u:
                    if isinstance(v, str):
                        if v == e_obj.v_alias:
                            return e_obj
                    else:
                        if v == e_obj.v:
                            return e_obj
        
        return None
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def find_fanin(self, u):
        """Returns incoming edge objects incident to u.

        Parameters
        ----------
        u : int or str
            Node identifier.
            If a string is passed, an alias-based lookup is performed.
    
        Returns
        -------
        List[Edge]
            The list of all incoming edges.
        """

        incident = []
        for e_obj in self.edges:
            if isinstance(u, str):
                if u == e_obj.v_alias:
                    incident.append(e_obj)
            else:
                if u == e_obj.v:
                    incident.append(e_obj)

        return incident
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def find_fanout(self, u):
        """Returns outgoing edge objects incident to u.

        Parameters
        ----------
        u : int or str
            Node identifier.
            If a string is passed, an alias-based lookup is performed.
    
        Returns
        -------
        List[Edge]
            The list of all outgoing edges.
        """

        incident = []
        for e_obj in self.edges:
            if isinstance(u, str):
                if u == e_obj.u_alias:
                    incident.append(e_obj)
            else:
                if u == e_obj.u:
                    incident.append(e_obj)

        return incident
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def alias_nodes(self, aliases):
        """Assigns aliases to (a subset of) nodes and relabels
        edges accordnigly.

        Parameters
        ----------
        aliases : Dict[int, str]
            Mapping of nodes to aliases.

        Returns
        -------
        None
        """

        for u in aliases:
            u_obj = self.find_node(u)
            if u_obj is not None:
                u_obj.set_alias(aliases[u])
            fanin = self.find_fanin(u)
            for i in fanin:
                i.set_head_alias(aliases[u])
            fanout = self.find_fanout(u)
            for o in fanout:
                o.set_tail_alias(aliases[u])
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def generate_pos_overlap_constraints(self):
        """Generates the exact position overlap removal constraints.
        
        Fixed nodes are assumed to be implicitly movable within their
        respective clusters. This is a necessary precondition for
        optimization by moving a small subset of nodes to succeed.
        If we want to reserve a position for a particular fixed node,
        we should do that by externally removing all conflicting cov-pairs.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        pos_dict = {}
        for u_obj in self.nodes:
            for pos in u_obj.positions:
                try:
                    pos_dict[pos].append(u_obj.u)
                except:
                    pos_dict.update({pos : [u_obj.u]})
        
        self.pos_overlap_constraints = []
        
        #String conversions:
        suffix = lambda pos : "_%d_%d_%d" % (pos.x, pos.y, pos.i)
        u_pos_var = lambda u, pos : "u_%d" % u + suffix(pos)

        for pos in pos_dict:
            if len(pos_dict[pos]) < 2:
                continue
            cst = ""
            for u in pos_dict[pos]:
                cst += u_pos_var(u, pos) + " + "
            cst = cst[:-len(" + ")] + " <= 1"
            self.pos_overlap_constraints.append(cst)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def generate_constraints(self):
        """Generates constraints for all objects.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        for e in self.edges:
            e.generate_constraints(self.fixed_nodes)
        for u in self.nodes:
            u.generate_constraints()
        for c in self.clusters:
            c.generate_constraints(self.rev_adj, self.fixed_nodes)

        self.generate_pos_overlap_constraints()

        self.timing.generate_constraints()

        if INCIDENCE_COUNTING:
            self.count_fanouts()
            self.count_fanins()
            self.infer_pattern_nodes()
            self.fanout_implications()
            self.fanin_implications()
            #Remove the position assignments for those
            #nodes for which the fanouts and fanins have
            #been removed after implication encoding.
            rem_list = []
            for u in self.pat_pos_assigns:
                if not u in self.fo_assigns and not u in self.fi_assigns:
                    rem_list.append(u)
            for u in rem_list:
                del self.pat_pos_assigns[u]
                for p in range(0, N):
                    try:
                        self.bin_vars.remove("u_%d_%d" % (u, p))
                    except:
                        pass

        #Add the "if no cov, no pos" constraints for each node.
        #For this we need a global view.

        #String conversions:

        pos_suffix = lambda pos : "%d_%d_%d" % (pos.x, pos.y, pos.i)
        pos_var = lambda u, pos : "u_%d_%s" % (u, pos_suffix(pos))
        edge_var = lambda e, cov : "e_%d_%d_%s_%s"\
                 % (e.u, e.v, pos_suffix(cov.u), pos_suffix(cov.v))

        pos_dict = {}
        for e in self.cov_map:
            u, v = e
            for cov in self.cov_map[e]:
                try:
                    pos_dict[pos_var(u, cov.u)].add(edge_var(e, cov))
                except:
                    pos_dict.update({pos_var(u, cov.u) : set([edge_var(e, cov)])})
                try:
                    pos_dict[pos_var(v, cov.v)].add(edge_var(e, cov))
                except:
                    pos_dict.update({pos_var(v, cov.v) : set([edge_var(e, cov)])})

        self.only_cov_pos = ""
        for pos in pos_dict:
            only_cov_pos = "- " + pos + " + "
            for e in pos_dict[pos]:
                only_cov_pos += e + " + "
            only_cov_pos = only_cov_pos[:-len(" + ")] + " >= 0\n"
            self.only_cov_pos += only_cov_pos
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

        txt = "Subject to\n"
        if FORCE_COVERAGE_INCREASE:
            txt += "\\* Minimum Coverage \\*\n"
            txt += self.least_coverage
        txt += "\***   Edges   ***\\\n"
        for e in self.edges:
            txt += e.print_constraints()
            if not LAZY_STRETCHES:
                txt += e.print_lazy_constraints() 
        txt += "\***   Nodes   ***\\\n"
        for u in self.nodes:
            txt += u.print_constraints()
            if not LAZY_STRETCHES:
                txt += u.print_lazy_constraints()
        txt += self.only_cov_pos
        txt += "\***   Clusters   ***\\\n"
        for c in self.clusters:
            txt += c.print_constraints()
        txt += "\***   Positions   ***\\\n"
        for p in self.pos_overlap_constraints:
            txt += p + "\n"
        txt += "\***   Timing   ***\\\n"
        if REMOVE_IRRELEVANT_EDGES:
            max_delays = {}
            for e in self.edges:
                #NOTE:We assume that the edges were appropriately aliased.
                u_str = e.u_alias
                if u_str.endswith(".d"):
                    u_str = u_str[:-1] + 'q'
                v_str = e.v_alias
                if v_str.endswith(".q"):
                    v_str = v_str[:-1] + 'd'
                e_str = (u_str, v_str)
                max_delays.update({e_str : self.maad[e_str]})
            self.timing.remove_irrelevant_edges(max_delays)
        txt += self.timing.print_constraints()

        if LAZY_STRETCHES:
            txt += "\nLAZY CONSTRAINTS\n"
            for e in self.edges:
                txt += e.print_lazy_constraints()
            for u in self.nodes:
                txt += u.print_lazy_constraints()

        if INCIDENCE_COUNTING:
            txt += "\***   Degree-matching   ***\\\n"
            for u in self.fo_assigns:
                txt += self.fo_assigns[u] + "\n"
            for u in self.fi_assigns:
                txt += self.fi_assigns[u] + "\n"
            for u in self.pat_pos_assigns:
                txt += self.pat_pos_assigns[u]
            for u in self.fanout_implications:
                txt += self.fanout_implications[u]
            for u in self.fanin_implications:
                txt += self.fanin_implications[u]

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

        txt = "Bounds\n"

        if FILTER_SLOW:
            txt += "\***   Edge Delays   ***\\\n"
            for e in self.edges:
                txt += e.print_td_bound()

        #TODO: See if a full reduction of the sum to achieve complete
        #      binarization of variables would be useful.
        txt += "\***   Input Indicators   ***\\\n"
        for c in self.clusters:
            txt += c.print_input_indicator_bounds()
        txt += "\***   Timing   ***\\\n"
        if TIGHT_TA_BOUNDS:
            max_delays = {}
            min_delays = {}
            for e in self.edges:
                #NOTE:We assume that the edges were appropriately aliased.
                u_str = e.u_alias
                if u_str.endswith(".d"):
                    u_str = u_str[:-1] + 'q'
                v_str = e.v_alias
                if v_str.endswith(".q"):
                    v_str = v_str[:-1] + 'd'
                e_str = (u_str, v_str)
                max_delays.update({e_str : e.td_bound_num[1]})
                min_delays.update({e_str : e.td_bound_num[0]})
            positive_slacks = self.timing.tighten_arrival_bounds(min_delays, max_delays)
            if not positive_slacks:
                return None

        txt += self.timing.print_bounds(self.forced_cpd_lb)

        if INCIDENCE_COUNTING:
            txt += "\***   Degree-bounds   ***\\\n"
            for u in self.fo_bounds:
                txt += self.fo_bounds[u] + "\n"
            for u in self.fi_bounds:
                txt += self.fi_bounds[u] + "\n" 

        return txt
    #-----------------------------------------------------------------------#
    
    #-----------------------------------------------------------------------#
    def print_variables(self):
        """Prints all variables.

        NOTE: Continuous variables need not be declared.

        Parameters
        ----------
        None
    
        Returns
        -------
        str
            Text of the variables.
        """

        for e in self.edges:
            self.bin_vars |= e.bin_vars
        for u in self.nodes:
            self.bin_vars |= u.bin_vars
        for c in self.clusters:
            self.bin_vars |= c.bin_vars

        for c in self.clusters:
            self.general_vars |= c.general_vars

        txt = "Binary\n"
        for v in sorted(self.bin_vars):
            txt += v + "\n"
        txt += "General\n"
        for v in sorted(self.general_vars):
            txt += v + "\n"

        return txt
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def print_ilp(self):
        """Prints the entire ILP formulation.

        Parameters
        ----------
        None

        Returns
        -------
        str
            Text of the ILP.
        """

        empty_coverage = False
        if OPT_MODE == "opt":
            txt = "Minimize Tamax\n"
        elif MAXIMIZE_COVERAGE_SIMILARITY:
            if not hasattr(self, "incumbent_coverage"):
                print("set_incumbent_coverage must be called.")
                raise ValueError
            txt = "Maximize %s \n" % self.incumbent_coverage
            if self.incumbent_coverage == "NOTHING":
                empty_coverage = True
        elif MINIMIZE_DISPLACEMENT:
            txt = "Minimize "
            for u in self.nodes:
                txt += u.displacement
            txt = txt[:-len(" + ")] + "\n"
        elif OPTIMIZER == "glpk":
            txt = "Minimize NOTHING\n"
        else:
            txt = "Minimize Tamax\n"
        txt += self.print_constraints()
        bounds = self.print_bounds()
        if bounds is None:
            return None

        txt += bounds
        if hasattr(self, "coverage_fixing"):
            txt += self.coverage_fixing
        if OPT_MODE == "feas" and (OPTIMIZER == "glpk" or empty_coverage):
            #NOTE: CPLEX can accept 0 directly in the objective function,
            #but GLPK can not.
            txt += "NOTHING = 0\n"
        txt += self.print_variables()
        txt += "End\n"
        
        return txt
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def solve_cplex(self, timeout):
        """Calls CPLEX to solve the ILP.

        Parameters
        ----------
        timeout : int
            Time allowed for solving, in the number of seconds.

        Returns
        -------
        None
        """

        with open("ilp.lp", "w") as prob:
            txt = self.print_ilp()
            if txt is None:
                return
            prob.write(txt)
   
        call = "cplex -c \"set timelimit %d\" " % timeout
        call += "\" set mip strategy search 1\" "
        call += "\" set mip strategy heuristicfreq -1\" "
        if OPT_MODE == "feas":
            call += "\"set emphasis mip 1\" "
            call += "\" set mip limits solutions 1\" "
        call += "\"read ilp.lp\" "
        if self.read_mst:
            call += "\"read ilp_save.mst\" "
        call += "\"optimize\" \"set logfile ilp.sol\" "
        call += "\"display solution variables -\" \"write ilp.mst\" \"quit\""
        call += " > cplex_call.log"
        os.system(call)
        os.system("mv ilp.mst ilp_save.mst")
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def solve_glpk(self, timeout):
        """Calls GLPK to solve the ILP.

        Parameters
        ----------
        timeout : int
            Time allowed for solving, in the number of seconds.

        Returns
        -------
        None
        """

        with open("ilp.lp", "w") as prob:
            txt = self.print_ilp()
            if txt is None:
                return
            prob.write(txt)

        call = "glpsol --tmlim %d --lp ilp.lp" % timeout
        call += " --write ilp.sol --wglp ilp.prob"
        call += " > glpsol.log"

        os.system(call)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def solve(self, timeout):
        """Calls the appropriate solver to sove the ilp.
        
        Parameters
        ----------
        timeout : int
            Time allowed for solving, in the number of seconds.

        Returns
        -------
        None
        """

        if OPTIMIZER == "glpk":
            self.solve_glpk(timeout)
        elif OPTIMIZER == "cplex":
            self.solve_cplex(timeout)
        else:
            print("Unknown optimizer.")
            exit(-1)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def parse_solution_cplex(self):
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

        try:
            with open("ilp.sol", "r") as sol:
                lines = sol.readlines()
        except:
            return {}

        val_dict = {}
        rd = False
        zero_others = False
        for line in lines:
            if "Variable Name" in line:
                rd = True
                continue
            if not rd:
                continue
            if "All other variables" in line and "are 0" in line:
                zero_others = True
                break
            var = line.split()[0]
            val = line.split()[1] 
            if var in self.bin_vars or var in self.general_vars:
                val = int(float(val))
                #NOTE: CPLEX dumps everything as float.
            else:
                val = float(val)
            val_dict.update({var : val})
        
        if zero_others:   
            for var in self.bin_vars:
                if not var in val_dict:
                    val_dict.update({var : 0})
            for var in self.general_vars:
                if not var in val_dict:
                    val_dict.update({var : 0})
            #NOTE: Maybe you want to add continuous variables as well.

        self.solution = val_dict

        return val_dict
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def parse_solution_glpk(self):
        """Parses the GLPK-produced solution, returning a dictionary of values
        of all variables.

        Parameters
        ----------
        None

        Returns
        -------
        Dict[str, float]
            Dictionary of variable values.
        """

        col_dict = {}
        try:
            with open("ilp.prob", "r") as prob:
                lines = prob.readlines()
        except:
            return {}

        for line in lines:
            words = line.split()
            if len(words) != 4 or words[1] != 'j':
                continue
            col_dict.update({int(words[2]) : words[3]})

        with open("ilp.sol", "r") as sol:
            lines = sol.readlines()

        val_dict = {}
        for line in lines:
            words = line.split()
            if len(words) != 3 or words[0] != 'j':
                continue
            var = col_dict[int(words[1])]
            if var in self.bin_vars or var in self.general_vars:
                val = int(words[2])
            else:
                val = float(words[2])
            val_dict.update({var : val})

        self.solution = val_dict

        return val_dict
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def parse_solution(self):
        """Parses the solution, returning a dictionary of values
        of all variables.

        Parameters
        ----------
        None

        Returns
        -------
        Dict[str, float]
            Dictionary of variable values.
        """

        if OPTIMIZER == "glpk":
            return self.parse_solution_glpk()
        elif OPTIMIZER == "cplex":
            return self.parse_solution_cplex()
        else:
            print("Unknown optimizer.")
            exit(-1)
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def get_arrival_times(self):
        """Returns all arrival times.

        Parameters
        ----------
        None
    
        Returns
        -------
        Dict[str, float]
            A dictionary of node arrival times.

        Raises
        ------
        AssertionError
            If there is a variable for a node with fixed arrival time.
            If any arrival time is negative.
            If arrival time of any node is smaller than the arrival
            times of its parents.
        """

        #Fetchers:
        fetch_node = lambda u : self.inv_node_counts[int(u.split("ta_")[1])]
        
        tarrs = {}
        for var in self.solution:
            try:
                u = fetch_node(var)
            except:
                continue
            val = self.solution[var]
            tarrs.update({u : val})

        if __debug__:
            for u in tarrs:
                u_tar = tarrs[u]
                for p in self.tg.pred[u]:
                    try:
                        p_tar = tarrs[p]
                    except:
                        p_tar = self.tarrs[p]
                    assert (p_tar <= u_tar),\
                           "(%s (%d), %s(%d)): Arrival times inverted"\
                            % (p, self.node_counts[p], u, self.node_counts[u])

        self.sol_tarrs = tarrs

        return tarrs
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def get_edge_delays(self, tolerance):
        """Creates a copy of the timing graph and annotates
        its edges with the delays from the ILP solution.

        Also performs a consistency check.

        Parameters
        ----------
        tolerance : float
            Floating point comparison tolerance.

        Returns
        -------
        None
        
        Raises
        ------
        AssertionError
            If tolerance is not a nonnegative float.
            If any edge delay is not equal to one of the covering /
            stretching options.
        """
        
        check_nonneg_float(tolerance)

        #Fetchers:
        fetch_u = lambda e : self.inv_node_counts[int(e.split("td_")[1]\
                             .split('_')[0])]
        fetch_v = lambda e : self.inv_node_counts[int(e.split("td_")[1]\
                             .split('_')[1])]

        tg = self.tg.copy()
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

            if __debug__:
                e = BasicEdge(self.node_counts[u], self.node_counts[v])
                possible = []
                try:
                    possible += [cov.td for cov in self.cov_map[e]]
                except:
                    pass
                try:
                    possible += [stretch.td for stretch in self.stretch_map[e]]
                except:
                    pass
                assert (possible and any(abs(td - val) < tolerance\
                                         for td in possible)),\
                       "Found impossible edge delay."

            tg[u][v]["td"] = val

        self.sol_tg = tg
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def check_solution_timing_consistency(self, tolerance):
        """Checks if the arrival times reported in the solution are
        consistent with the chosen edge delays.

        Parameters
        ----------
        tolerance : float
            Floating point comparison tolerance.

        Returns
        -------
        None
        
        Raises
        ------
        AssertionError
            If tolerance is not a nonnegative float.
        LookupError
            If an arrival time can not be found in the solution or among
            the fixed arrival times.
        ValueError
            If any arrival time is not consistent with the computed one.
            If there is any inconsistency among edge delays (i.e., if
            the current positions of the endpoint nodes dictate a delay
            different from the one found in the solution).
        """
       
        check_nonneg_float(tolerance)

        tg = self.sol_tg
        for v in nx.topological_sort(self.sol_tg):
            tg.node[v]["ta"] = max([tg.node[u]["ta"] + tg[u][v]["td"]\
                                    for u in tg.pred[v]] + [0])
            ta = -1
            try:
                ta = self.sol_tarrs[v]
            except:
                try:
                    ta = self.tarrs[v]
                except:
                    print("Arrival time not found in either the solution or"\
                          + " among the fixed arrival times: %s (%d)"\
                          % (v, self.node_counts[v]))
                    raise LookupError
            #Arrival times in the solution can be larger than the real ones,
            #but they can not be smaller. We already checked that the solution
            #arrival times are mutually consistent (i.e., that each node has
            #solution arrival time at least as large as all of its parents).
            if tg.node[v]["ta"] > ta + tolerance:
                print(v + "(%d) Arrival times inconsistent." % self.node_counts[v])
                print("Calculated: %.4f, In solution: %.4f" % (tg.node[v]["ta"], ta))
                raise ValueError

        #Check edge delay consistency.
        for u, v in tg.edges():
            found = False
            e = BasicEdge(self.node_counts[u], self.node_counts[v])
            if u in self.pos_dict:
                if v in self.pos_dict:
                    for cov in self.cov_map.get(e, []):
                        if cov.u == self.pos_dict[u]\
                           and cov.v == self.pos_dict[v]:
                            found = True
                            if abs(cov.td - tg[u][v]["td"]) > tolerance:
                                print("Edge delay inconsistency.")
                                raise ValueError
                            break
                    if not found:
                        u_cls = node_pos_to_node_cls(self.pos_dict[u])
                        v_cls = node_pos_to_node_cls(self.pos_dict[v])
                        for stretch in self.stretch_map.get(e, []):
                            if stretch.u == u_cls and stretch.v == v_cls:
                                found = True
                                if abs(stretch.td - tg[u][v]["td"]) > tolerance:
                                    print("Edge delay inconsistency.")
                                    raise ValueError
                                break
                else:
                    u_cls = node_pos_to_node_cls(self.pos_dict[u])
                    try:
                        v_cls = self.cls_dict[v]
                    except:
                        v_cls = tg.node[v]["coords"]
                    for stretch in self.stretch_map.get(e, []):
                        if stretch.u == u_cls and stretch.v == v_cls:
                            found = True
                            if abs(stretch.td - tg[u][v]["td"]) > tolerance:
                                print("Edge delay inconsistency.")
                                raise ValueError
                            break
            else:
                try:
                    u_cls = self.cls_dict[u]
                except:
                    u_cls = tg.node[u]["coords"]
                try:
                    v_cls = node_pos_to_node_cls(self.pos_dict[v])
                except:
                    try:
                        v_cls = self.cls_dict[v]
                    except:
                        v_cls = tg.node[v]["coords"]
                for stretch in self.stretch_map.get(e, []):
                    if stretch.u == u_cls and stretch.v == v_cls:
                        found = True
                        if abs(stretch.td - tg[u][v]["td"]) > tolerance:
                            print("Edge delay inconsistency.")
                            raise ValueError
                        break
            if not found:
                if abs(tg[u][v]["td"] - self.tg[u][v]["td"]) > tolerance:
                    print("Edge delay inconsistency: %s -> %s (%d -> %d)"\
                          % (u, v, self.node_counts[u], self.node_counts[v]))
                    print("Warning: Not covered by any cov or stretch.")
                    print("Expected: %.4f, In solution: %.4f"\
                          % (self.tg[u][v]["td"], tg[u][v]["td"]))
                    raise ValueError
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def get_positions(self):
        """Returns the positions of all movable nodes

        Also performs a consistency check.

        Parameters
        ----------
        None

        Returns
        -------
        Dict[str, NodePos]
            Exact positions for those nodes on which they are enforced.
        Dict[str, NodeCls]
            Clusters, for those nodes without an exact position.

        Raises
        ------
        ValueError
            If there are position overlaps.
            If a node is assigned multiple exact positions.
            If a node is assigned zero or > 1 clusters.
            If there is a mismatch between the position and the cluster.
            If the given positions do not imply the assigned edge delay.
        """

        #Fetchers:
        #.......................................................................#
        def fetch_pos(self, var):
            """Retrurns a node and its position implied by the variable.

            Parameters
            ----------
            var : str
                Variable name.
            
            Returns
            -------
            str
                Node identifier.
            NodePos
                Position.
            """

            if not var.startswith("u_"):
                return None, None
            suffix = var.split("u_", 1)[1]
            terms = suffix.split('_')
            if len(terms) != 4:
                return None, None

            u = self.inv_node_counts[int(terms[0])]
            pos = NodePos(int(terms[1]), int(terms[2]), int(terms[3]))
            
            return u, pos
        #.......................................................................#

        #.......................................................................#
        def fetch_cls(self, var):
            """Retrurns a node and its cluster implied by the variable.

            Parameters
            ----------
            var : str
                Variable name.
            
            Returns
            -------
            str
                Node identifier.
            NodeCls
                Cluster.
            """

            if not var.startswith("u_"):
                return None, None
            suffix = var.split("u_", 1)[1]
            terms = suffix.split('_')
            if len(terms) != 3:
                return None, None

            u = self.inv_node_counts[int(terms[0])]
            cls = NodeCls(int(terms[1]), int(terms[2]))
            
            return u, cls
        #.......................................................................#

        pos_dict = {}
        cls_dict = {}
        for var in self.solution:
            if self.solution[var] == 0:
                continue
            u, pos = fetch_pos(self, var)
            if pos is None:
                u, cls = fetch_cls(self, var)
                if cls is None:
                    continue
                if u in cls_dict:
                    print(u + "(%s): Reassigning cluster." % var)
                    raise ValueError
                cls_dict.update({u : cls})
                if u.endswith(".d"):
                    cls_dict.update({u[:-1] + 'q' : cls})
                elif u.endswith(".q"):
                    cls_dict.update({u[:-1] + 'd' : cls})
                if u in pos_dict and node_pos_to_node_cls(pos_dict[u]) != cls:
                    print(u + ": Position and cluster mismatch.")
                    raise ValueError
            else:
                if u in pos_dict:
                    print(u + ": Reassigning position.")
                    raise ValueError
                pos_dict.update({u : pos})
                if u.endswith(".d"):
                    pos_dict.update({u[:-1] + 'q' : pos})
                elif u.endswith(".q"):
                    pos_dict.update({u[:-1] + 'd' : pos})
                if cls_dict.get(u, None) not in (None, node_pos_to_node_cls(pos)):
                    print(u + ": Position and cluster mismatch.")
                    raise ValueError

        #Strip the clusters for those nodes that have an exact position.
        for u in pos_dict:
            try:
                del cls_dict[u]
            except:
                pass

        present = list(pos_dict.keys()) + list(cls_dict.keys())
        if not all(u_obj.alias in present for u_obj in self.nodes):
            print("Not all nodes assigned a position or a cluster.")
            raise ValueError

        self.pos_dict = pos_dict
        self.cls_dict = cls_dict

        ff_strip = lambda u : u[:-len(".d")] if u.endswith((".d", ".q"))\
                              else u

        pos_usage = {p : set([(ff_strip(u), self.node_counts[u])\
                         for u in pos_dict if pos_dict[u] == p])\
                     for p in pos_dict.values()}

        for p in pos_usage:
            if len(pos_usage[p]) > 1:
                print("Position", p, "overused.")
                print pos_usage[p]
                raise ValueError

        return pos_dict, cls_dict
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def reconstruct_coverage(self):
        """Returns a list of edges that were successfully covered,
        labeled with the original node labels.

        Parameters
        ----------
        None

        Returns
        -------
        Dict[Tuple[str], PhysDirCon]
            A dictionary with all covered connections listing the physical
            connections that were used for covering.
        """


        #Fetchers:
        
        u_num = lambda cov : int(cov.split('_')[2])
        v_num = lambda cov : int(cov.split('_')[3])
        u_pos = lambda cov : NodePos(int(cov.split('_')[3]),\
                                     int(cov.split('_')[4]),\
                                     int(cov.split('_')[5]))
        v_pos = lambda cov : NodePos(int(cov.split('_')[6]),\
                                     int(cov.split('_')[7]),\
                                     int(cov.split('_')[8]))
        cov_con = lambda cov : PhysDirCon(u_pos(cov), v_pos(cov))

        cov_vars = [var for var in self.bin_vars if var.startswith("cov_ind")\
                    and self.solution[var]]
        e_vars = [var for var in self.bin_vars if var.startswith("e_")\
                  and len(var.split('_')) == 9 and self.solution[var]]
        cov_dict = {}
        for cov in cov_vars:
            u = self.inv_node_counts[u_num(cov)]
            if u.endswith(".d"):
                u = u[:-1] + 'q'
            v = self.inv_node_counts[v_num(cov)]
            if v.endswith(".q"):
                v = v[:-1] + 'd'
            con = None
            for e in e_vars:
                if e[1:].startswith(cov[len("cov_ind"):]):
                    con = cov_con(e)
                    break
            cov_dict.update({(u, v) : con})

        return cov_dict
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def fix_coverage(self, coverage):
        """Fixes the appropriate cov-indicators if a coverage list is
        supplied, or even exact physical connection mapping, if a
        dictionary holding it is supplied.

        Parameters
        ----------
        coverage : List[Tuple[str]] or Dict[Tuple[str], PhysDirCon]
            Coverage to be enforced.

        Returns
        -------
        None
        """

        phys_suffix = lambda con : "_%d_%d_%d_%d_%d_%d"\
                    % (con.u.x, con.u.y, con.u.i, con.v.x, con.v.y, con.v.i)
        e_prefix = lambda u, v : "_%d_%d" % (self.node_counts[u], self.node_counts[v])

        self.coverage_fixing = ""
        if isinstance(coverage, list) or isinstance(coverage, set):
            for e in coverage:
                self.coverage_fixing += "cov_ind" + e_prefix(e[0], e[1]) + " = 1\n"
        elif isinstance(coverage, dict):
            for e in coverage:
                self.coverage_fixing += 'e' + e_prefix(e[0], e[1])\
                                      + phys_suffix(coverage[e]) + " = 1\n"
    #-----------------------------------------------------------------------#


    #-----------------------------------------------------------------------#
    def check_cov_consistency(self, pattern):
        """Checks the cov selection for consistency with node positions
        and then checks if this that selection really forms a subgraph of
        the periodic pattern graph.

        Parameters
        ----------
        pattern : List[FreeDirCon]

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If any inconsistencies occur.
        """

        fanin_dict = {}
        fanout_dict = {}
        for e in self.cov_map:
            u_str = self.inv_node_counts[e.u]
            if u_str.endswith(".d"):
                u_str = u_str[:-1] + 'q'
            v_str = self.inv_node_counts[e.v]
            if v_str.endswith(".q"):
                v_str = v_str[:-1] + 'd'
            try:
                fanin_dict[v_str].append(u_str)
                fanout_dict[u_str].append(v_str)
            except:
                fanin_dict.update({v_str : [u_str]})
                fanout_dict.update({u_str : [v_str]})

        #String conversions:

        cov_ind = lambda u, v : "cov_ind_%d_%d"\
                  % (self.node_counts[u], self.node_counts[v])

        fanout_cov_dict = {}
        for u in fanout_dict:
            for v in fanout_dict[u]:
                if self.solution[cov_ind(u, v)]:
                    u_pos = self.pos_dict.get(u, None)
                    if u_pos is None:
                        print("Node", u, "(%d)" % self.node_counts[u],\
                              "not assigned a fixed position.")
                        raise ValueError
                    v_pos = self.pos_dict.get(v, None)
                    if v_pos is None:
                        print("Node", v, "(%d)" % self.node_counts[v],\
                              "not assigned a fixed position.")
                        raise ValueError
                    x_offset = v_pos.x - u_pos.x
                    y_offset = v_pos.y - u_pos.y
                    u_ind = u_pos.i
                    v_ind = v_pos.i
                    found = False
                    for e in pattern:
                        if e.u == u_ind and e.x_offset == x_offset\
                           and e.v == v_ind and e.y_offset == y_offset:
                            found = True
                            break
                    if not found:
                        print("Covering of %s -> %s (%d -> %d) not possible"\
                              % (u, v, self.node_counts[u], self.node_counts[v])\
                              + " with u_pos = (%d, %d, %d) and v_pos = (%d, %d, %d)"\
                              % (u_pos.x, u_pos.y, u_pos.i, v_pos.x, v_pos.y, v_pos.i))
                        raise ValueError
                    try:
                        fanout_cov_dict[u].append(v)
                    except:
                        fanout_cov_dict.update({u : [v]})
        fanin_cov_dict = {}
        for v in fanin_dict:
            for u in fanin_dict[v]:
                if self.solution[cov_ind(u, v)]:
                    u_pos = self.pos_dict.get(u, None)
                    if u_pos is None:
                        print("Node", u, "(%d)" % self.node_counts[u],\
                              "not assigned a fixed position.")
                        raise ValueError
                    v_pos = self.pos_dict.get(v, None)
                    if v_pos is None:
                        print("Node", v, "(%d)" % self.node_counts[v],\
                              "not assigned a fixed position.")
                        raise ValueError
                    x_offset = v_pos.x - u_pos.x
                    y_offset = v_pos.y - u_pos.y
                    u_ind = u_pos.i
                    v_ind = v_pos.i
                    found = False
                    for e in pattern:
                        if e.u == u_ind and e.x_offset == x_offset\
                           and e.v == v_ind and e.y_offset == y_offset:
                            found = True
                            break
                    if not found:
                        print("Covering of %s -> %s (%d -> %d) not possible"\
                              % (u, v, self.node_counts[u], self.node_counts[v])\
                              + " with u_pos = (%d, %d, %d) and v_pos = (%d, %d, %d)"\
                              % (u_pos.x, u_pos.y, u_pos.i, v_pos.x, v_pos.y, v_pos.i))
                        raise ValueError
                    try:
                        fanin_cov_dict[v].append(u)
                    except:
                        fanin_cov_dict.update({v : [u]})

        pattern_fanin_dict = {}
        pattern_fanout_dict = {}
        for e in pattern:
            try:
                pattern_fanout_dict[e.u] += 1
            except:
                pattern_fanout_dict.update({e.u : 1})
            if not e.v in fanout_dict:
                fanout_dict.update({e.v : 0})
            try:
                pattern_fanin_dict[e.v] += 1
            except:
                pattern_fanin_dict.update({e.v : 1})
            if not e.u in fanin_dict:
                fanin_dict.update({e.u : 0})
 
        for u in fanin_cov_dict:
            indeg = len(fanin_cov_dict[u])
            u_pos = self.pos_dict[u].i
            pos_indeg = pattern_fanin_dict[u_pos]
            if indeg > pos_indeg:
                print("In-degree of node", u, "at position", u_pos, "too high.")
                print("Required: %d, Available: %d" % indeg, pos_indeg)
                raise ValueError
     
        for u in fanout_cov_dict:
            outdeg = len(fanout_cov_dict[u])
            u_pos = self.pos_dict[u].i
            pos_outdeg = pattern_fanout_dict[u_pos]
            if outdeg > pos_outdeg:
                print("Out-degree of node", u, "at position", u_pos, "too high.")
                print("Required: %d, Available: %d" % outdeg, pos_outdeg)
                raise ValueError
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def generate_branching_order(self):
        """Generates a branching priority order.
        """

        order = {"tau" : (2, "UP"), "e" : (1, "UP")}
        #NOTE: Try first branching on tau, which fixes the cov/uncov indicators,
        #that in turn fix the fanout/fanin counters and many position variables,
        #then on e-variables (cov-pairs) which fix node positions, and in turn
        #clusters. Everything else is not really relevant. In both cases, branch
        #first up, as that will zero most of other variables, because of the 
        #at-most-one and exactly-one constraints.

        #ORD format: 1) First line = NAME
        #            2) Last line = ENDATA
        #            3) all other lines = space + UP (DN) + space + name + space + priority

        ordered_vars = {t : [] for t in order}
        for var in self.bin_vars:
            for t in order:
                if var.startswith(t):
                    ordered_vars[t].append(var)
        txt = "NAME\n"
        for t in sorted(order, key = lambda t : order[t][0], reverse = True):
            for var in sorted(ordered_vars[t]):
                line = " %s %s %d\n" % (order[t][1], var, order[t][0])
                txt += line
        txt += "ENDATA"

        with open("ilp.ord", "w") as dumpf:
            dumpf.write(txt)
    #-----------------------------------------------------------------------#
##########################################################################
