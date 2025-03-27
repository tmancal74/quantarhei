# -*- coding: utf-8 -*-

from ..symbolic.cumulant import Uop
from ..symbolic.cumulant import UopEater
from ..symbolic.cumulant import transform_to_einsum_expr

class DSFeynmanDiagram():
    """ Double-sided Feynman diagrams


    """
    
    def __init__(self, ptype="non-defined"):

        self.states = {}
        self.pointer = 0
        self._pic_rep = "\n"
        self.finished = False

        self.type = ptype
        
        self.states[self.pointer] = ["g", "g"]
        self.add_states_line()
        
        self.dimensions = {"t1":0, "t2":-1, "t3":1}
        
        self.diag_name = None


    def add_states_line(self):
        self._pic_rep = "      | "+self.states[self.pointer][0]+ \
                        "      "+self.states[self.pointer][1]+" |\n"+self._pic_rep


    def add_arrow(self, side, dir, to="a"):
        """Adds an arrow in a given direction and transits to a new state

        """
        self.pointer += 1
        if side == "left":
            self.states[self.pointer] = [to, self.states[self.pointer-1][1]]
            if dir == "--->":
                self._pic_rep = "  --->|----------|\n"+self._pic_rep
            elif dir == "<---":
                self._pic_rep = "  <---|----------|\n"+self._pic_rep
                
        elif side == "right":
            self.states[self.pointer] = [self.states[self.pointer-1][0],to]
            if dir == "<---":
                self._pic_rep = "      |----------|<---\n"+self._pic_rep
            elif dir == "--->":
                self._pic_rep = "      |----------|--->\n"+self._pic_rep
                
        self.add_states_line()


    def finish(self, end="g"):
        """Finishes the diagram with a given state

        """
        self.add_arrow("left", "<---", end)
        self._pic_rep = "\n"+self._pic_rep
        self.count = self.pointer - 1
        self.finished = True


    def _check_finished(self):
        """Checkes if the diagram is finished

        """
        if not self.finished:
            raise Exception("An unfinisged diagram: finish() method has to be called first")


    def get_time_dimensions(self):
        """ Returns a dictionary describing the required representation 
        of the response as a matrix of times. 
        
        """
        return self.dimensions

    
    def get_phase_factor(self, dimensions=None):
        """Returns the phase factor for the present diagram, with reshaped time symbols if provided."""
    
        fact = "-"
        Nst = len(self.states)
    
        def reshape_time_symbol(tname):
            """Return reshaped time variable like t1[:,None] or t2[None,None]"""
            if dimensions is None:
                return tname
    
            axis = dimensions.get(tname, None)
            if axis == 0:
                return f"{tname}[:,None]"
            elif axis == 1:
                return f"{tname}[None,:]"
            elif axis == -1:
                return f"{tname}[None,None]"
            else:
                raise ValueError(f"Unsupported or missing axis for time variable '{tname}' in dimensions.")
    
        for sec in range(1, Nst - 1):
            st = self.states[sec]
            fact += "1j*(En["+st[0]+"]-En["+st[1]+"])"
            
            if dimensions is None:
                tm = "t"+str(sec)
                fact += "*"+tm+" "
            else:
                # Get time name by index
                tm = list(dimensions.keys())[sec - 1]
                tm_reshaped = reshape_time_symbol(tm)
                fact += "*"+tm_reshaped+" "
    
            if sec < Nst - 2:
                fact += "-"
    
        return fact        


    def evolution_operators(self, operators=False):
        """Returns evolution operators corresponding to the diagram
        
        It can return a string representation or a list of operator objects.
        
        Parameters
        ----------
        
        operators : bool
            If operators is False, string representation will be returned. 
            A list of Uop objects will be returned if operators is True.
        
        """

        self._check_finished()
            
        kk = 0
        Uops = ""
        Uops_list = []
 
        for key in self.states:
            if kk > 0 and kk < self.count+1:
                sts = self.states[key]
                rightstate = sts[1]
                
                times = dict()
                if kk == 1:  
                    times["t1"] = 1
                elif kk == 2:
                    times["t2"] = 1 
                elif kk == 3:
                    times["t3"] = 1 
                else:
                    raise Exception("Unknown time for evolution operator")
                    
                if rightstate == "g":
                    Uops += "Ugd(t"+str(kk)+")"
                    if operators:
                        Uops_list.append(Uop(state="g", times=times, dagger=True))
                    
                else:
                    Uops += "Ued("+rightstate+",t"+str(kk)+")"
                    if operators:
                        Uops_list.append(Uop(state=rightstate, times=times, dagger=True))
                        
                if kk < self.count:
                    Uops += "*"
            kk += 1

        kk = 0
        rUop = ""
        rUops_list = []
        for key in self.states:
             if kk > 0 and kk < self.count+1:
                sts = self.states[key]
                leftstate = sts[0]
                
                times = dict()
                if kk == 1:  
                    times["t1"] = 1
                elif kk == 2:
                    times["t2"] = 1 
                elif kk == 3:
                    times["t3"] = 1 
                else:
                    raise Exception("Unknown time for evolution operator")                
                
                if leftstate == "g":
                    rUop = "Ug(t"+str(kk)+")"+rUop
                    if operators:
                        list_st = [Uop(state="g", times=times, dagger=False)]
                        list_st.extend(rUops_list)
                        rUops_list = list_st
                else:
                    rUop = "Ue("+leftstate+",t"+str(kk)+")"+rUop
                    if operators:
                        list_st = [Uop(state=leftstate, times=times, dagger=False)]
                        list_st.extend(rUops_list)
                        rUops_list = list_st
                        
                if kk < self.count:
                    rUop = "*"+rUop
             kk += 1

        Uops += "*"+rUop
        Uops_list.extend(rUops_list)
        
        if operators:
            return Uops_list
        
        return Uops
    
    
    def coherence_GF(self):
        """Returns coherence Green's function product for this diagram 
        
        """
        
        evs = self.evolution_operators(operators=True)
        eater = UopEater()
        out_list = eater.eat(evs)
        return eater.spit_coherence_GF()
    
    
    def get_cumulant_expression(self, verbose=False):
        """Returns the cumulant evaluation of the diagram
        
        """

        codes = []

        code_import = """
import sympy as sp
import quantarhei as qr

from quantarhei.symbolic.cumulant import Ugde, Uedg, Uged, Uegd
from quantarhei.symbolic.cumulant import evaluate_cumulant

"""

        codes.append(code_import)

        code_symbols = """
a = sp.Symbol("a")
b = sp.Symbol("b")

t1 = sp.Symbol("t1")
t2 = sp.Symbol("t2")
t3 = sp.Symbol("t3")

"""
        outs = self.coherence_GF()
    
        code1 = code_symbols+"A_cum = "+outs
        codes.append(code1)

        code2 = "\n"+\
            "evc = evaluate_cumulant(A_cum, positive_times=[t1, t2, t3])\n"
#        code2 = "\n"+\
#            "evc = evaluate_cumulant(A_cum)\n"
       
        codes.append(code2)        
       
        local_vars = {}
        for cod in codes:
            if verbose:
                print(cod)
            compiled_code = compile(cod, "<string>","exec")
            exec(compiled_code, {}, local_vars)
        
        return local_vars["evc"]   
    
    
    def __str__(self):
        self._check_finished()
        return self._pic_rep

    
    def report(self):
        self._check_finished()
        print("\nDiagram of",self.type,"type\n")
        print(self)
        print("Light interaction count:", self.count)   


    def get_vectorized_code(self, function=True, participation_matrix=None):
        """ Return the code that evaluates the response function
        
        """
        
        dims = self.dimensions
        phfac = self.get_phase_factor(dimensions=dims)
        #print("\n ... phase factor:", phfac)

        coh = self.coherence_GF()
        #print("\n ... in coherence Green's functions:", coh)

        cme = self.get_cumulant_expression()
        #print("\nCumulants:\n", cme)

        # TASK: the expression with participation_matrix must include a correct energy phase factor
        #       This factor must be constructed based on the cumulant expression if possible. If not,
        #       we need to submit some information on it separately.
        if participation_matrix is not None:
            prt = participation_matrix
            out = transform_to_einsum_expr(cme,
                                participation_matrix=participation_matrix,
                                dimensions=dims)
        else:
            out = transform_to_einsum_expr(cme,
                                       participation_matrix=None,
                                       dimensions=dims)
            
        out_str = str(out)+" "+phfac
        
        if function:
            
            fcode = _format_code(3, out_str)

            #prt = "MM"
            fstr = "\ndef "+self.diag_name+"(t2, t1, t3, system):"
            
            fstr += \
'''
    """ Returns a matrix of the respose function values for given t1 and t3 
    
    Parameters:
    -----------
    
    t1 : numpy.array
        Array of t1 times (must be the same as the t1 axis of the gg object)
        
    t2 : float
        Value of the t2 (waiting) time of the response
        
    t3 : numpy.array
        Array of t3 times (must be the same as the t3 axis of the gg object)
        
    system : aggregate or molecule class
        An object storing all information about the system including 
        the values of the line shape functions.
    
    
    """'''
            
            fstr += "\n    import numpy as np"
            
            fstr += "\n"
            # FIXME: Molecule and Aggregate has to have something like this
            
            fstr += "\n    gg = system.get_lineshape_functions()"
            # FIXME: FunctionStorage should be able to return this
            if participation_matrix is not None:
                fstr += "\n    # Mx = system.get_participation()"
                fstr += "\n    "+prt+" = system.get_weighted_participation()"
            # FIXME: make sure something like this can be obtained from aggregate and molecule
            fstr += "\n    En = system.get_eigenstate_energies()"
            fstr += "\n    g = 1  # ground state index"
            fstr += "\n    gg.create_data(reset={'t2':t2})"
            fstr += "\n"
            fstr += "\n    Ne = En.shape[0]"
            fstr += "\n    ret = numpy.zeros((len(t1),len(t3)), dtype=COMPLEX)"
            
            # FIXME: Here we have to allow diffent number of loops
            fstr += "\n    for a in range(Ne):"
            fstr += "\n        for b in range(Ne):"
            fstr += "\n"
            fstr += "\n            ret += \\\n"
            fstr += fcode
            fstr += "\n"
            fstr += "\n    return ret"
            
            return fstr
            
        return out_str
    

def _format_code(N, code_string):
    """
    Formats a long Python expression string into a properly indented, multiline expression
    wrapped inside numpy.exp(...), with line breaks only at top-level '+' and '-' operators.

    Parameters:
    - N: number of 4-space indents
    - code_string: expression to wrap and format

    Returns:
    - A formatted string ready for Python code
    """
    indent = " " * (4 * N)
    inner_indent = indent + "    "

    def smart_split(expr):
        parts = []
        current = ""
        depth = 0
        in_quote = False
        quote_char = ""

        i = 0
        while i < len(expr):
            char = expr[i]

            # Handle quotes to avoid splitting inside strings
            if char in {"'", '"'}:
                if in_quote and char == quote_char:
                    in_quote = False
                    quote_char = ""
                elif not in_quote:
                    in_quote = True
                    quote_char = char

            if not in_quote:
                if char in "([{":
                    depth += 1
                elif char in ")]}":
                    depth -= 1

                # Split at top-level + or - signs (skip unary - at start)
                if depth == 0 and i > 0 and expr[i-1] != "e" and char in "+-":
                    parts.append(current.strip())
                    current = char
                    i += 1
                    continue

            current += char
            i += 1

        if current.strip():
            parts.append(current.strip())
        return parts

    # Remove line breaks from input
    code_string = code_string.replace("\n", "")
    parts = smart_split(code_string)

    if not parts:
        return indent + "np.exp(0)"

    # Construct formatted output
    formatted = indent + "np.exp(\n"
    for part in parts:
        formatted += inner_indent + part + "\n"
    formatted += indent + ")"
    return formatted


class R1g_Diagram(DSFeynmanDiagram):
    """R1g diagram

    Diagram of R1g type
    
    
          | g      g |
      <---|----------|
          | a      g |
          |----------|--->
          | a      b |
          |----------|<---
          | a      g |
      --->|----------|
          | g      g |

    """
    
    def __init__(self, states=["a","b"]):
        super().__init__(ptype="R1g")
        self.add_arrow("left", "--->", to=states[0])
        self.add_arrow("right", "<---", to=states[1])
        self.add_arrow("right", "--->", "g" )
        self.finish()
        
        self.diag_name="R1g"



class R2g_Diagram(DSFeynmanDiagram):
    """R2g diagram

    Diagram of R2g type
    
    
          | g      g |
      <---|----------|
          | b      g |
          |----------|--->
          | b      a |
      --->|----------|
          | g      a |
          |----------|<---
          | g      g |

    """

    def __init__(self, states=["a","b"]):
        super().__init__(ptype="R2g")
        self.add_arrow("right", "<---", to=states[0])
        self.add_arrow("left", "--->", to=states[1])
        self.add_arrow("right", "--->", "g" )
        self.finish()

        self.diag_name="R2g"


class R3g_Diagram(DSFeynmanDiagram):
    """R3g diagram

    Diagram of R3g type
    
    
          | g      g |
      <---|----------|
          | b      g |
      --->|----------|
          | g      g |
          |----------|--->
          | g      a |
          |----------|<---
          | g      g |

    """

    def __init__(self, states=["a","b"]):
        super().__init__(ptype="R3g")
        self.add_arrow("right", "<---", to=states[0])
        self.add_arrow("right", "--->", to="g")
        self.add_arrow("left", "--->", to=states[1])
        self.finish()
    

class R4g_Diagram(DSFeynmanDiagram):
    """R4g diagram

    Diagram of R4g type
    
    
          | g      g |
      <---|----------|
          | b      g |
      --->|----------|
          | g      g |
      <---|----------|
          | a      g |
      --->|----------|
          | g      g |

    """

    def __init__(self, states=["a","b"]):
        super().__init__(ptype="R4g")
        self.add_arrow("left", "--->", to=states[0])
        self.add_arrow("left", "<---", to="g")
        self.add_arrow("left", "--->", to=states[1])
        self.finish()


class R1f_Diagram(DSFeynmanDiagram):
    """R1f diagram

    Diagram of R1f type
    
    
          | b      b |
      <---|----------|
          | f      b |
      --->|----------|
          | a      b |
          |----------|<---
          | a      g |
      --->|----------|
          | g      g |

    """

    def __init__(self, states=["a","b","f"]):
        super().__init__(ptype="R1f")
        self.add_arrow("left", "--->", to=states[0])
        self.add_arrow("right","<---", to=states[1])
        self.add_arrow("left", "--->", to=states[2])
        self.finish(end="b")
        

class R2f_Diagram(DSFeynmanDiagram):
    """R2f diagram

    Diagram of R2f type
    
          | a      a |
      <---|----------|
          | f      a |
      --->|----------|
          | b      a |
      --->|----------|
          | g      a |
          |----------|<---
          | g      g |
          
    """

    def __init__(self, states=["a","b","f"]):
        super().__init__(ptype="R2f")
        self.add_arrow("right", "<---", to=states[0])
        self.add_arrow("left","--->", to=states[1])
        self.add_arrow("left", "--->", to=states[2])
        self.finish(end="a")
