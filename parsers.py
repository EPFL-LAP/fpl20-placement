"""Module contains all VPR output parsers needed by the ILP placer.
"""

import networkx as nx

from feeder_types import Ble
from ilp_placer_types import NodeCls, FreeDirCon, FreeProgCon

from config_vars import N

##########################################################################
class PlacementTimingGraph(object):
    """Models the postplacement timing graph.

    Parameters
    ----------
    tg_filename : str
        Name of the timing graph file.
    net_filename : str
        Name of the packed netlist file.

    Arguments
    ---------
    *Parameters
    tg : nx.DiGraph
        Basic timing graph, simply as parsed from the input file.

    Methods
    -------
    parse_net(mem : bool, mult : bool, io : bool)
        Parses the packed netlist to obtain the order of appearance of
        signal identifiers.
    """

    #-----------------------------------------------------------------------#
    def __init__(self, tg_filename, net_filename):
        """Constructor of the ParseTimingGraph class.
        """
        
        self.tg_filename = tg_filename
        self.net_filename = net_filename
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def parse_net(self, mem, mult, io):
        """Parses the packed netlist to obtain the order of appearance of
        signal identifiers.

        Nodes of the VPR timing graph appear in the same order as in
        the packed netlist. To be able to annotate them with signal
        identifiers, we must obtain that order.

        Parameters
        ----------
        mem : bool
            Include memory pins.
        mult : bool
            Include multiplier pins.
        io : bool
            Include I/O pins.

        Returns
        -------
        None
        """

        #Initialization:
        sigs = {"seq" : {}, "comb" : {}, "in" : {}, "out" : {}}

        #FFs and memory ports are sequential, LUTs, multiplier ports, and
        #repeaters are combinational, and I/Os have a separate timing
        #graph node category.
        
        memory_blocks = set()
        port_inits = ["<inputs>", "<outputs>", "<clocks>"]
        port_ends = ["</inputs>", "</outputs>", "</clocks>"]
        state = "idle"
        open_blocks = 0
        block_name_stack = []
        blif_name = ""

        #Fetchers:
        fetch_attr = lambda l, a : l.split(a + "=\"", 1)[1].split('"', 1)[0]
        fetch_ind = lambda s : int(s.split('[', 1)[1].split(']', 1)[0])
        fetch_name = lambda l : fetch_attr(l, "name")
        fetch_mode = lambda l : fetch_attr(l, "mode")
        fetch_inst = lambda l : fetch_attr(l, "instance")
        fetch_assignments = lambda l : l.split('>', 1)[1].split('<', 1)[0].split()
        fetch_last = lambda s : s.rsplit("___", 1)[-1]

        #.......................................................................#
        def merge_assignments(lines):
            """Merges together assignments broken accross multiple lines
            into one signle line.
            """
            
            merged = []
            stick = False
            nl = ""
            for line in lines:
                nl += line
                if "<port " in line and not "</port>" in line:
                    stick = True
                if "</port>" in line:
                    stick = False
                if not stick:
                    merged.append(nl)
                    nl = ""

            return merged
        #.......................................................................#
        with open(self.net_filename, "r") as net:
            lines = merge_assignments(net.readlines())

        for lcnt, line in enumerate(lines[1:], start = 1):
            if state == "idle":
                if "<block" in line:
                    state = "port_init"
                    block_name = fetch_inst(line)
                    try:
                        block_mode = fetch_mode(line) 
                    except:
                        block_mode = ""
                    block_ind = fetch_ind(block_name)
                    sigs["seq"].update({block_ind : []})
                    sigs["comb"].update({block_ind : []})
                    if io and block_name.startswith("io["):
                        io_name = fetch_name(line)
                        io_mode = fetch_mode(line)
                        if  io_mode == "inpad":
                            sigs["in"].update({block_ind : [io_name]})
                        elif io_mode == "outpad":
                            sigs["out"].update({block_ind : [io_name]})
                    block_name_stack.append(block_name)
                    open_blocks += 1
            elif state == "port_init":
                if line.strip() in port_inits:
                    state = "port"
                elif "<block" in line:
                    if not "/>" in line:
                        blif_name = fetch_name(line)
                        block_name += "___" + fetch_inst(line)
                        block_name_stack.append(block_name)
                        open_blocks += 1
                    try:
                        block_mode = fetch_mode(line)
                    except:
                        block_mode = ""
                elif "</block" in line:
                    open_blocks -= 1
                    block_name_stack.pop(-1)
                    if open_blocks == 0:
                        state = "idle"
                        continue
                    block_name = block_name_stack[-1]
            elif state == "port":
                if line.strip() in port_ends:
                    state = "port_init"
                    continue
                port_name = fetch_name(line)
                added = 0
                pins = []
                for i, a in enumerate(fetch_assignments(line)):
                    if a != "open":
                        pins.append(block_name + '.' + port_name + '[' + str(i) + ']')
                        added = 1
                if added:
                    #FIXME: Make this more architecture independent.
                    if fetch_last(pins[-1]) == "lut[0].out[0]":
                        sigs["comb"][block_ind].append(blif_name)
                    elif 'Q' in fetch_name(line):
                        sigs["seq"][block_ind].append(blif_name)
                    elif "name=\"out\">repeater_" in line:
                        sigs["comb"][block_ind].append(blif_name)
                    elif mult and block_mode == ""\
                         and any(fetch_last(block_name).startswith(m)\
                                 for m in ["mult_36x36[", "mult_18x_18[", "mult9x9["]):
                        if "out" in fetch_name(line):
                            for w in fetch_assignments(line):
                                if w != "open":
                                    sigs["comb"][block_ind].append(w)
                    elif mem and "memory[" in block_name and not "___" in block_name:
                        memory_blocks.add(block_ind)
                        if any(n in fetch_name(line) for n in ["addr", "data", "we"]):
                            sigs["seq"][block_ind] += pins
                    elif mem and "memory_slice[" in block_name:
                        memory_blocks.add(block_ind)
                        if "out" in fetch_name(line):
                            for w in fetch_assignments(line):
                                if w != "open":
                                    sigs["seq"][block_ind].append(w)
    
        self.net_sigs = sigs
        self.memory_blocks = memory_blocks
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def parse_tg(self):
        """Parses the timing graph and annotates the nodes
        according to the signal identifier sequence produced
        by parse_net.

        NOTE: We assume that the node coordinates are embedded in the
              timing graph, as in our modification of VTR-7.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Raises
        ------
        AssertionError
            If any of the lists in self.net_sigs is not empty at the end
            of the parsing process. This indicates that either a bug is
            present in the parser, or, more likely, that the timing graph
            file does not correspond to the given packed netlist file.
        """
        
        #Initialization:        
        nodes = []
        edges = []
        blif_nodes = {}
        coords = {}
        reading = 0

        #Fetchers:
        fetch_source = lambda l : int(l.split()[0].split('(')[0])
        fetch_x = lambda l : int(l.split()[0].split('(')[1].split(',')[0])
        fetch_y = lambda l : int(l.split()[0].split('(')[1].split(',')[1][:-1])
        fetch_coords = lambda l : NodeCls(fetch_x(l), fetch_y(l))
        fetch_block_ind = lambda l : int(line.split()[3])
        next_sig = lambda b_ind, t : self.net_sigs[t][b_ind].pop(0)
        next_comb = lambda b_ind : next_sig(b_ind, "comb")
        next_seq_d = lambda b_ind : next_sig(b_ind, "seq") + ".d"
        next_seq_q = lambda b_ind : next_sig(b_ind, "seq") + ".q"
        next_in = lambda b_ind : next_sig(b_ind, "in")
        next_out = lambda b_ind : next_sig(b_ind, "out")

        with open(self.tg_filename, "r") as tg:
            lines = tg.readlines()

        for lcnt, line in enumerate(lines):
            if line.isspace():
                continue
            if "TN" in line:
                u = fetch_source(line)
                nodes.append(u)
                coords.update({u : fetch_coords(line)})
                block_ind = fetch_block_ind(line)
                if "CB_IPIN" in line and block_ind in self.memory_blocks:
                    #Tie the memory input drivers to CB_IPINs.
                    if self.net_sigs["seq"].get(block_ind, [""])[0]\
                                               .startswith("memory"):
                        #We have to skip the clk pin as we did not store it.
                        u_blif = next_seq_d(block_ind)
                        assert (u_blif not in blif_nodes), "Updating nodes."
                        blif_nodes.update({u_blif : u})
                if "OUTPAD_SINK" in line:
                    try:
                        u_blif = next_out(block_ind)
                        assert (u_blif not in blif_nodes)
                        blif_nodes.update({u_blif : u})
                    except AssertionError:
                        print("Updating nodes.")
                    except:
                        pass
                elif "INPAD_SOURCE" in line:
                    try:
                        u_blif = next_in(block_ind)
                        assert (u_blif not in blif_nodes)
                        blif_nodes.update({u_blif : u})
                    except AssertionError:
                        print("Updating nodes.")
                    except:
                        pass
                if any(s in line for s in ["PRIMITIVE_OP", "CONSTANT_GEN"]):
                    try:
                        u_blif = next_comb(block_ind)
                        assert (u_blif not in blif_nodes)
                        blif_nodes.update({u_blif : u})
                    except AssertionError:
                        print("Updating nodes.")
                    except:
                        pass
                elif "FF_SINK" in line:
                    if block_ind in self.memory_blocks:
                        #For memories, we tie the driver pins to CB_IPINs
                        #instead of the sinks. This changes nothing in the
                        #analysis, as the periphery will take care of the rest,
                        #but saves us the trouble of distributing the driver
                        #onto all appropriate sinks (VPR creates a unique sink
                        #for each output and each input that influences it
                        #(i.e., each pair of such)).
                        continue
                    try:
                        ff_d = next_seq_d(block_ind)
                        assert (ff_d not in blif_nodes)
                        blif_nodes.update({ff_d : u})
                        #We must save the FF for the Q pin as well.
                        self.net_sigs["seq"][block_ind].insert(0, ff_d[:-len(".d")])
                    except AssertionError:
                        print("Updating nodes.")
                    except:
                        pass
                elif "FF_SOURCE" in line:
                    skip_ind = 0
                    if block_ind in self.memory_blocks:
                        while self.net_sigs["seq"].get(block_ind, [""])[skip_ind]\
                                                      .startswith("memory"):
                            #Skip all the sink-only ports.
                            skip_ind += 1
                    try:
                        u_blif = next_seq_q(block_ind)
                        assert (u_blif not in blif_nodes)
                        blif_nodes.update({u_blif : u})
                    except AssertionError:
                        print("Updating nodes.")
                    except:
                        pass
                if any(s in line for s in ["FF_SINK", "FF_CLOCK", "OUTPAD_SINK"]):
                    #These pins have no outgoing edges
                    #(in fact, the SINK/SOURCE pair opens the FFs)
                    reading = 0
                    continue

                v = int(line.split()[-2])
                td = float(line.split()[-1])
                if abs(td) > 1e1:
                    #This is a constant generator. Set its delay to zero.
                    td = 0
                reading = 1
                edges.append((u, v, {"td" : td}))
            elif "num_tnode_levels" in line:
                break
            elif reading:
                v = int(line.split()[0])
                td = float(line.split()[-1])
                edges.append((u, v, {"td" : td}))
   
        assert (all(all(len(self.net_sigs[cat][b]) == 0 for b in self.net_sigs[cat])\
                for cat in self.net_sigs)), "Timing graph and .net do not match."
 
        tg = nx.DiGraph()
        tg.add_edges_from(edges) 

        self.tg = tg
        self.blif_nodes = blif_nodes
        self.coords = coords
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def relabel_nodes(self):
        """Relabels the timing graph nodes and coordinate keys to
        self.blif_nodes.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        
        relabeling_dict = {self.blif_nodes[u] : u for u in self.blif_nodes}
        self.tg = nx.relabel_nodes(self.tg, relabeling_dict)

        coords = {}
        for u in self.coords:
            try:
                coords.update({relabeling_dict[u] : self.coords[u]})
            except:
                coords.update({u : self.coords[u]})

        self.coords = coords
    #-----------------------------------------------------------------------#

    #-----------------------------------------------------------------------#
    def merge_coords(self):
        """Merges coordinates to node attributes.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        for u in self.coords:
            self.tg.node[u]["coords"] = self.coords[u]
    #-----------------------------------------------------------------------#
##########################################################################

##########################################################################
def parse_bles(net_filename, cluster = False, ble_keyword = "fle"):
    """Returns the list of all BLEs in the packed netlist.

    Parameters
    ----------
    net_filename : str
        Name of the packed netlist file.
    cluster : Optional[bool]
        Determines if the BLEs should be returned in a flat
        list or clustered into CLB lists.
    
    Returns
    -------
    List[Ble]
        The list of BLEs.
    or List[List[Ble]]
        Depending on the value of >>cluster<<
    """

    #Fetchers:
    fetch_attr = lambda l, a : l.split(a + "=\"", 1)[1].split('"', 1)[0]
    fetch_name = lambda l : fetch_attr(l, "name")
    fetch_inst = lambda l : fetch_attr(l, "instance")
    fetch_mode = lambda l : fetch_attr(l, "mode")
    is_lut = lambda l : fetch_inst(l) == "lut[0]"\
                   or fetch_inst(l) == "lut6[0]" and fetch_mode(l) == "wire" 
    is_ff = lambda l : fetch_inst(l) == "ff[0]"
    is_open_ble = lambda l : fetch_name(l) == "open"\
                  and fetch_inst(l).startswith(ble_keyword)

    with open(net_filename, "r") as net:
        lines = net.readlines()

    open_blocks = 0
    prev_open_blocks = 0
    bles = []
    clustered_bles = []
    skip = True
    for line in lines:
        if "<block" in line:
            block_name = fetch_inst(line)
            if open_blocks == 1:
                if block_name.startswith("clb["):
                    skip = False
                    clustered_bles.append([])
                else:
                    skip = True
                    continue
            blif_node = fetch_name(line)
            lut_stored = is_lut(line)
            ff_stored = is_ff(line)
            if "/>" in line:
                if lut_stored:
                    clustered_bles[-1].append([blif_node])
                    bles.append([blif_node])
                elif ff_stored:
                    #NOTE: We assume that in a BLE with an FF
                    #there is always a LUT in wire mode as well.
                    clustered_bles[-1][-1].append(blif_node)
                    bles[-1].append(blif_node)
                elif is_open_ble(line):
                    clustered_bles[-1].append(["open", "open"])
            else:
                prev_open_blocks = open_blocks
                open_blocks += 1
        elif not skip and "</block>" in line:
            if prev_open_blocks < open_blocks:
                if lut_stored:
                    clustered_bles[-1].append([blif_node])
                    bles.append([blif_node])
                elif ff_stored:
                    clustered_bles[-1][-1].append(blif_node)
                    bles[-1].append(blif_node)
            prev_open_blocks = open_blocks
            open_blocks -= 1

    erase = lambda s : "" if s == "open" else s 
    if cluster:
        for i in range(0, len(clustered_bles)):
           for cnt, ble in enumerate(clustered_bles[i]):
               clustered_bles[i][cnt] = Ble(erase(ble[0]), erase(ble[1]))
        return clustered_bles
    else:
        for cnt, ble in enumerate(bles):
            bles[cnt] = Ble(erase(ble[0]), erase(ble[1]))
        return bles
##########################################################################

##########################################################################
def parse_grid_types(grid_types_filename):
    """Reads the FPGA grid cell type map, as output by our
    modification of VTR-7

    Parameters
    ----------
    grid_types_filename : str
        Name of the grid types map file.

    Returns
    -------
    Dict[NodeCls, str]
        A dictionary of types indexed by cell coordinates.
    """ 

    with open(grid_types_filename, "r") as grid_types:
        lines = grid_types.readlines()

    grid_types = {}
    for line in lines:
        words = line.split()
        x = int(words[0])
        y = int(words[1])
        tp = words[2]
        grid_types.update({NodeCls(x, y) : tp})

    return grid_types
##########################################################################

##########################################################################
def parse_pattern(pattern_filename):
    """Parses a pattern.

    Parameters
    ----------
    pattern_filename : str
        Name of the pattern file.

    Returns
    -------
    List[FreeDirCon]
        Pattern represented by a list of static graph edges.
    """

    with open(pattern_filename, "r") as pattern:
        pattern_text = pattern.read()

    pattern = []
    state = "idle"
    for c in pattern_text:
        if c.isspace():
            continue
        if c == '#':
            state = "comment"
        if state == "comment":
            if c == '}':
                state = "idle"
            else:
                continue
        elif state == "idle" and c == '{':
            pattern.append([])
            state = "tail"
            tail = ""
        elif state == "tail":
            if c == ',':
                state = "head"
                head = ""
                pattern[-1].append(tail)
            else:
                tail += c
        elif state == "head":
            if c == ',':
                state = "offset_x"
                offset_x = ""
                pattern[-1].append(head)
            else:
                head += c
        elif state == "offset_x":
            if c == '(':
                pattern[-1].append([])
            elif c == ',':
                state = "offset_y"
                offset_y = ""
                offset_x = int(offset_x)
                pattern[-1][-1].append(offset_x)
            else:
                offset_x += c
        elif state == "offset_y":
            if c == ')':
                state = "id"
                edge_id = ""
                offset_y = int(offset_y)
                pattern[-1][-1].append(offset_y)
            else:
                offset_y += c
        elif state == "id":
            if c == ',':
                continue
            if c == '}':
                state = "idle"
                pattern[-1].append(edge_id)
            else:
                edge_id += c
        else:
            print "Unknown state."
            exit(-1)

    #Convert to the required format.
    #We assume that in the input file nodes are specified as
    # LUT[u].PORT[pin_high:pin_low]
    fetch_node = lambda u : int(u.split('.', 1)[0].split('[')[1][:-1])

    pattern_conv = []
    for e in pattern:
        pattern_conv.append(FreeDirCon(fetch_node(e[0]), fetch_node(e[1]),\
                                       e[2][0], e[2][1]))

    return pattern_conv
##########################################################################

##########################################################################
def parse_vpr_delay_table(delay_table_filename):
    """Parses the VPR place-time intercluster delay lookup table.

    We expect that the table is dumped as in our modification of VTR-7.

    Parameters
    ----------
    delay_table_filename : str
        Name of the delay table file.

    Returns
    -------
    Dict[FreeDirCon, float]
        Dictionary of connection delays,
        indexed by cluster location pairs.
    """

    #.......................................................................#
    def read_one_table(start_line, end_line):
        """Reads a single of the multiple tables in the file.

        Parameters
        ----------
        start_line : str
            Line marking the start of the table.
        end_line : str
            Line marking the end of the table.
        
        Returns
        -------
            Dict[FreeDirCon, float]
            Dictionary of connection delays,
            indexed by cluster location pairs.
        """
         
        dt = []
        rd = False
        for line in lines:
            if line.strip() == start_line:
                rd = True
            elif line.strip() == end_line:
                break
            if rd:
                words = line.split()
                try:
                    td = float(words[0])
                except:
                    continue
                dt.append([])
                for w in words:
                    dt[-1].append(float(w))
        dt.reverse()
        dt_dict = {}
        for y, row in enumerate(dt):
            for x, td in enumerate(row):
                dt_dict.update({FreeProgCon(x, y) : td})
    
        return dt_dict
    #.......................................................................#

   
    with open(delay_table_filename, "r") as dt:
        lines = dt.readlines()

    start_line = "printing delta_clb_to_clb_sym"
    end_line = "printing delta_clb_to_clb_qpp"
    clb_clb = read_one_table(start_line, end_line)

    start_line = "printing delta_io_to_clb"
    end_line = "printing delta_clb_to_io"
    io_clb = read_one_table(start_line, end_line)

    start_line = "printing delta_clb_to_io"
    end_line = "printing delta_io_to_io"
    clb_io = read_one_table(start_line, end_line)

    start_line = "printing delta_io_to_io"
    end_line = "FILE_END"
    io_io = read_one_table(start_line, end_line)

    return {"clb" : {"clb" : clb_clb, "io" : clb_io},\
            "io" : {"clb" : io_clb, "io" : io_io}}
##########################################################################

##########################################################################
def parse_pattern_delays(pattern, delays_filename):
    """Parses the pattern connection delays.

    The input format is: "CLB-offset fractional-CLB-offset delay",
    where we assume that the LUTs are stacked vertically in the CLB,
    as in the Stratix series of FPGAs.

    Parameters
    ----------
    pattern : List[FreeDirCon]
        Pattern represented by a list of static graph edges.
    delays_filename : str
        Name of the delays file.
    
    Returns
    -------
    Dict[FreeDirCon, float]
        A dictionary of delays of the pattern connections.
    """

    #........................................................................#
    def compute_offsets(e):
        """Computes the integral and fractional CLB offsets for an edge.

        Parameters
        ----------
        e : FreeDirCon
            Edge for which offsets are being computed.
        
        Returns
        -------
        Tuple[int]
            A pair of offsets.
        """

        y_offset = e.y_offset    
        if y_offset == 0:
            return abs(e.x_offset), abs(e.v - e.u)

        if y_offset < 0:
            y_offset += 1
            fract_offset = N + e.u - e.v
            if fract_offset >= N:
                y_offset -= 1
                fract_offset -= N
        else:
            y_offset -= 1
            fract_offset = N + e.v - e.u  
            if fract_offset >= N:
                y_offset += 1
                fract_offset -= cluster_size
    
        return abs(e.x_offset) + abs(y_offset), fract_offset
    #........................................................................#


    #Fetchers:
    fetch_int  = lambda l : int(l.split()[0])
    fetch_fract = lambda l : int(l.split()[1])
    fetch_delay = lambda l : float(l.split()[2])

    with open(delays_filename, "r") as df:
        lines = df.readlines()

    delays = {}
    for line in lines:
        delays.update({(fetch_int(line), fetch_fract(line))\
                       : fetch_delay(line)})

    delays_dict = {}
    for e in pattern:
        int_offset, fract_offset = compute_offsets(e)
        delays_dict.update({e : delays[(int_offset, fract_offset)]})

    return delays_dict
##########################################################################
